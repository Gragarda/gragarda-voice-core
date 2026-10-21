#pragma once
#include <atomic>
#include <cstddef>
#include <memory>
#include <vector>

struct DenoiseState;
struct DFState;
struct kiss_fft_state;

// Coeur DSP independant d'OBS.
// trim -> HPF -> debruitage (RNNoise ou DeepFilterNet 3) -> de-reverb -> EQ soustractif -> compresseur -> de-esser
// -> EQ additif -> leveler (AGC lent, cible LUFS selon destination) -> limiteur.
// Reglages pilotes par 4 macros 0..1 ; 0.5 partout = reglages de reference valides a l'ecoute.

struct Biquad {
    float b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0;
    enum Type { HighPass, Peak, LowShelf, HighShelf, BandPass };
    void design(Type t, double sr, double f, double q, double gainDb);
};

class VoiceChain {
public:
    static constexpr int kMaxChannels = 8;
    static constexpr int kNumEq = 5;
    static constexpr int kDnFrame = 480; // trame RNNoise (10 ms a 48 kHz)

    enum Engine { Fast = 0, Max = 1 }; // Rapide = RNNoise, Qualite max = DeepFilterNet 3
    enum Destination { Stream = 0, Podcast = 1, VoiceOver = 2 };

    struct Macros {
        float clean = 0.5f, warmth = 0.5f, clarity = 0.5f, punch = 0.5f; // 0..1
        Destination dest = Stream;
        // Correction du micro (posee par la calibration, invisible) : ramene le micro vers une voix neutre, dB
        float corrLow = 0, corrPresence = 0, corrAir = 0;
        // De-reverb (pose par la calibration) : RT60 estime de la piece en secondes, 0 = hors service
        float rt60 = 0;
    };
    static float targetLufs(Destination d) { return d == Podcast ? -16.0f : d == VoiceOver ? -18.0f : -14.0f; }

    // Alloue et charge les modeles : hors thread audio. dfnModel = chemin du .tar.gz DeepFilterNet
    // (nullptr ou introuvable -> mode Max indisponible, retombe sur un simple retard).
    void prepare(double sampleRate, int numChannels, const char* dfnModel = nullptr);
    bool maxAvailable() const { return (bool)dfn_[0]; }
    void setGainDb(float db);
    // Avant/Apres : sort le signal brut, retarde de la meme latence et ramene au volume du traite
    // (comparaison honnete : sans saut de synchro ni "plus fort = meilleur").
    void setBypass(bool b) { bypass_ = b; }
    void setDenoise(bool b) { denoise_ = b; } // la latence ne change pas
    void setEngine(Engine e) { engine_ = e; }  // change la latence (cf. latencySamples)
    // Thread UI : applique au debut du prochain bloc audio.
    // ponytail: pending_ copie sans verrou, un bloc peut voir un melange ancien/nouveau ; sans consequence audible
    void setMacros(const Macros& m) { pending_ = m; dirty_ = true; }

    // Audio planaire : channels[c][i], float -1..1. Taille de bloc quelconque.
    void process(float* const* channels, int numChannels, size_t numFrames);

    // Latence FIXE en echantillons (a afficher pour compenser la camera)
    // FIFO 1 trame + retard interne du moteur (RNNoise 1 trame, DeepFilterNet 3 trames) + de-reverb 1 trame + limiteur
    int latencySamples() const { return latencySamples(engine_); }
    int latencySamples(Engine e) const { return (2 + (e == Max ? 3 : 1)) * kDnFrame + lookahead_; }
    double sampleRate() const { return sampleRate_; }

private:
    double sampleRate_ = 48000.0;
    int numChannels_ = 2;
    float targetGain_ = 1.0f;
    float currentGain_ = 1.0f;
    float smoothCoef_ = 0.0f; // lissage pour eviter les "zips" quand on tourne le potard
    bool bypass_ = false;
    bool denoise_ = true;
    Engine engine_ = Fast;
    Macros pending_;
    std::atomic<bool> dirty_{true};
    void applyMacros(const Macros& m);
    float dnWet_ = 1.0f; // Nettoyage < 0.5 : on remelange du signal non debruite
    int engineDelayFrames() const { return engine_ == Max ? 3 : 1; }

    // Debruitage : FIFO d'une trame (480) puis moteur ; hors service = meme retard via historique
    static constexpr int kHist = 4;
    struct DnFree { void operator()(DenoiseState* s) const; };
    struct DfFree { void operator()(DFState* s) const; };
    std::unique_ptr<DenoiseState, DnFree> dn_[kMaxChannels];
    std::unique_ptr<DFState, DfFree> dfn_[kMaxChannels];
    float dnIn_[kMaxChannels][kDnFrame] = {}, dnOut_[kMaxChannels][kDnFrame] = {};
    float dnHist_[kMaxChannels][kHist][kDnFrame] = {};
    int dnPos_ = 0, dnHistIdx_ = 0;
    void denoiseFrame(int nc);

    // De-reverb : STFT 960 (hop 480, fenetre racine de Hann), soustraction spectrale de la reverberation
    // tardive (Lebart 2001). Retard fixe d'une trame, actif ou non.
    static constexpr int kDrBins = kDnFrame + 1, kDrHist = 5; // 5 trames = 50 ms : debut de la queue tardive
    struct FftFree { void operator()(kiss_fft_state* s) const; };
    std::unique_ptr<kiss_fft_state, FftFree> fft_;
    float drWin_[2 * kDnFrame] = {};
    float drPrev_[kMaxChannels][kDnFrame] = {}, drOla_[kMaxChannels][kDnFrame] = {};
    std::vector<float> drPow_; // [kMaxChannels][kDrHist][kDrBins] puissance lissee des trames passees, alloue dans prepare
    int drIdx_ = 0;
    float drDecay_ = 0; // attenuation de la queue apres kDrHist trames ; 0 = hors service
    void dereverbFrame(int nc);

    // Filtres : HPF + EQ (meme coefficients pour tous les canaux, etat par canal)
    Biquad hpf_, eq_[kNumEq], dsSplit_, dsDetect_, corr_[3];
    float hpfZ_[kMaxChannels][2] = {}, eqZ_[kNumEq][kMaxChannels][2] = {}, corrZ_[3][kMaxChannels][2] = {};
    float dsSplitZ_[kMaxChannels][2] = {}, dsDetectZ_[kMaxChannels][2] = {};

    // Compresseur (stereo lie)
    float compEnvDb_ = -120.0f, compAtt_ = 0, compRel_ = 0;
    float compThreshDb_ = -20.0f, compRatio_ = 3.0f, compMakeup_ = 1.0f;
    // De-esser
    float dsEnv_ = 0, dsAtt_ = 0, dsRel_ = 0, dsThreshDb_ = -28.0f;
    // Leveler : sonie BS.1770 simplifiee (filtre K), court terme 3 s, gain lent
    Biquad kw1_, kw2_;
    float kwZ_[kMaxChannels][2][2] = {};
    float lvMom_ = 0, lvShort_ = 0, lvMomC_ = 0, lvShortC_ = 0;
    float lvGainDb_ = 0, lvTarget_ = -14.0f, lvStep_ = 0;
    // Avant/Apres : signal brut retarde + gain d'egalisation de volume
    static constexpr int kDryLen = 4096; // > latence max (2472)
    std::vector<float> dry_; // [kMaxChannels][kDryLen], alloue dans prepare
    int dryPos_ = 0;
    float abWet_ = 0, abDry_ = 0, abGain_ = 1.0f, abC_ = 0, abSmooth_ = 0;
    // Limiteur a anticipation
    int lookahead_ = 0;
    std::vector<float> delay_;   // [kMaxChannels][lookahead_]
    std::vector<float> reqGain_; // gain requis sur la fenetre d'anticipation
    int pos_ = 0;
    float limGain_ = 1.0f, limRel_ = 0;
};

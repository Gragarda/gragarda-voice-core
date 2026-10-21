#pragma once
#include <memory>
#include <vector>

struct kiss_fft_state;

// De-reverb commun au temps reel (VoiceChain) et au hors temps reel (gv_offline).
// STFT 960 (hop 480, fenetre racine de Hann), soustraction spectrale de la reverberation tardive (Lebart 2001) :
// tardive = puissance lissee d'il y a 50 ms, attenuee selon le RT60. Une instance par canal.
class Dereverb {
public:
    static constexpr int kFrame = 480, kBins = kFrame + 1, kHist = 5; // 5 trames = 50 ms : debut de la queue tardive

    void prepare(); // alloue : hors thread audio
    // RT60 de la piece en secondes ; <= 0.05 = hors service (simple retard)
    void setRt60(float s);
    // Une trame de kFrame echantillons, en place ; sortie retardee d'une trame, active ou non
    void process(float* frame);

private:
    struct FftFree { void operator()(kiss_fft_state* s) const; };
    std::unique_ptr<kiss_fft_state, FftFree> fft_;
    // alloues dans prepare (VoiceChain en garde 8 : pas sur la pile)
    std::vector<float> win_, prev_, ola_; // [2*kFrame], [kFrame], [kFrame]
    std::vector<float> pow_;              // [kHist][kBins] puissance lissee des trames passees
    int idx_ = 0;
    float decay_ = 0; // attenuation de la queue apres kHist trames ; 0 = hors service
};

// RT60 aveugle (s) d'un enregistrement de parole mono, 0 = non mesurable. floorDb = fond de bruit (dB de puissance).
// Une piece ne s'eteint pas plus vite que son RT60, la parole si : les chutes d'energie les plus rapides
// (fenetres de 80 ms) donnent la decroissance de la piece. Sous-evalue (~0.45 s mesure pour 0.8 s reel).
float estimateRt60(const std::vector<float>& x, double sampleRate, float floorDb);
// Au-dela : de-reverb active. 0.6 s estime = tres reverberant
constexpr float kRt60On = 0.3f, kRt60Red = 0.6f;

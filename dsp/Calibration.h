#pragma once
#include <atomic>
#include <cstddef>
#include <vector>

// Calibration : capture du signal brut (avant la chaine) puis analyse -> reglages recommandes.
// Independant d'OBS, comme VoiceChain.

// Tampon de capture mono. push() depuis le thread audio : sans verrou ni allocation.
class Capture {
public:
    void prepare(size_t maxSamples) { buf_.assign(maxSamples, 0.0f); }
    // Thread UI. Le thread audio ecrit ensuite jusqu'a 'samples' echantillons.
    void start(size_t samples) { size_.store(0); want_ = samples < buf_.size() ? samples : buf_.size(); active_.store(true); }
    void stop() { active_.store(false); }
    bool done() const { return size_.load() >= want_; }
    float progress() const { return want_ ? (float)size_.load() / (float)want_ : 0.0f; }
    // Thread audio : downmix mono
    void push(const float* const* ch, int nc, size_t n);
    // Thread UI, une fois done() : copie de ce qui a ete capture
    std::vector<float> take() const { return std::vector<float>(buf_.begin(), buf_.begin() + size_.load()); }

private:
    std::vector<float> buf_;
    std::atomic<size_t> size_{0};
    size_t want_ = 0;
    std::atomic<bool> active_{false};
};

enum class Light { Green, Orange, Red };

struct CalibrationReport {
    // Mesures (signal brut, avant Gain)
    float noiseLufs = -120, speechLufs = -120, snrDb = 0;
    float peakDb = -120, clipPct = 0;
    float lowDb = 0, presenceDb = 0, sibilanceDb = 0; // ecarts au profil moyen de voix, dB
    float sibPeakDb = 0;                              // pics de sifflantes au-dessus de la moyenne, dB
    float spreadDb = 0;                               // dynamique : ecart p95-p10 du court terme
    float rt60 = 0;                                   // reverberation estimee (aveugle), s ; 0 = non mesurable
    // Voyants
    Light noise = Light::Green, level = Light::Green, sibilance = Light::Green, clipping = Light::Green;
    Light reverb = Light::Green;
    // Reglages recommandes
    float gainDb = 0;
    int clean = 50, warmth = 50, clarity = 50, punch = 50;
    float corrLow = 0, corrPresence = 0, corrAir = 0; // correction du micro, dB
    float dereverbRt60 = 0;                           // -> Macros::rt60 ; 0 = de-reverb hors service
    bool valid = false; // false : pas assez de parole detectee
};

// noise : ~5 s de silence ; speech : ~15 s de lecture. Mono, sampleRate quelconque (48 kHz attendu).
CalibrationReport analyzeCalibration(const std::vector<float>& noise, const std::vector<float>& speech,
                                     double sampleRate);

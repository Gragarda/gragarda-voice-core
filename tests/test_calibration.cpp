// Auto-test de l'analyse de calibration, sans OBS.
#undef NDEBUG
#include "Calibration.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

static unsigned gSeed = 11;
static float noise1() { gSeed = gSeed * 1664525u + 1013904223u; return (gSeed >> 8) / 8388608.f - 1.f; }

// voix synthetique : harmoniques 1/h, f0 variable, syllabes 3 Hz ; 'tilt' en dB/octave sur les harmoniques
static std::vector<float> voice(float amp, double secs, float noiseAmp, float tilt = 0)
{
    const size_t n = (size_t)(secs * 48000);
    std::vector<float> x(n);
    for (size_t i = 0; i < n; ++i) {
        double t = i / 48000.0, f0 = 140 + 30 * std::sin(2 * 3.14159265 * 0.7 * t);
        double env = std::fmax(0, std::sin(2 * 3.14159265 * 3 * t)), v = 0;
        for (int h = 1; h < 60 && f0 * h < 20000; ++h)
            v += std::sin(2 * 3.14159265 * f0 * h * t + h) / h * std::pow(10.0, tilt * std::log2((double)h) / 20);
        x[i] = (float)(amp * env * v) + noiseAmp * noise1();
    }
    return x;
}

static std::vector<float> hiss(float amp, double secs)
{
    std::vector<float> x((size_t)(secs * 48000));
    for (auto& v : x) v = amp * noise1();
    return x;
}

static void show(const char* name, const CalibrationReport& r)
{
    std::printf("%-10s bruit %.1f parole %.1f SNR %.1f crete %.1f | grave %+.1f pres %+.1f sib %+.1f pics %.1f dyn %.1f"
                " -> gain %+.1f N%d C%d Cl%d P%d\n",
                name, r.noiseLufs, r.speechLufs, r.snrDb, r.peakDb, r.lowDb, r.presenceDb, r.sibilanceDb,
                r.sibPeakDb, r.spreadDb, r.gainDb, r.clean, r.warmth, r.clarity, r.punch);
}

int main()
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);

    // Silence total : rapport invalide, pas de crash
    {
        std::vector<float> z(48000 * 5, 0.f), s(48000 * 15, 0.f);
        assert(!analyzeCalibration(z, s, 48000).valid);
    }

    // Micro propre, niveau correct
    auto clean = analyzeCalibration(hiss(0.0003f, 5), voice(0.1f, 15, 0.0003f), 48000);
    show("propre", clean);
    assert(clean.valid && clean.noise == Light::Green && clean.clipping == Light::Green);
    assert(clean.snrDb > 40);

    // Micro bruyant : Nettoyage monte, voyant bruit rouge
    auto noisy = analyzeCalibration(hiss(0.02f, 5), voice(0.1f, 15, 0.02f), 48000);
    show("bruyant", noisy);
    assert(noisy.noise == Light::Red && noisy.clean > clean.clean);

    // Micro faible : Gain positif ; micro fort : Gain negatif
    auto quiet = analyzeCalibration(hiss(0.0001f, 5), voice(0.02f, 15, 0.0001f), 48000);
    auto loud = analyzeCalibration(hiss(0.001f, 5), voice(0.5f, 15, 0.001f), 48000);
    show("faible", quiet);
    show("fort", loud);
    assert(quiet.gainDb > 10 && loud.gainDb < 0);
    assert(std::fabs((quiet.speechLufs + quiet.gainDb) - (loud.speechLufs + loud.gainDb)) < 0.5f);

    // Saturation : voyant rouge
    auto sat = voice(0.1f, 15, 0.0003f);
    for (auto& v : sat) v = std::fmax(-1.f, std::fmin(1.f, v * 20));
    auto clip = analyzeCalibration(hiss(0.0003f, 5), sat, 48000);
    show("sature", clip);
    assert(clip.clipping == Light::Red);

    // Voix sombre vs brillante : Clarte et Chaleur compensent dans le bon sens
    auto dark = analyzeCalibration(hiss(0.0003f, 5), voice(0.1f, 15, 0.0003f, -4), 48000);
    auto bright = analyzeCalibration(hiss(0.0003f, 5), voice(0.1f, 15, 0.0003f, +3), 48000);
    show("sombre", dark);
    show("brillante", bright);
    assert(dark.corrPresence > bright.corrPresence && dark.corrAir > bright.corrAir); // sombre eclairci, brillant adouci
    assert(dark.clarity == 50 && bright.clarity == 50); // macros timbre = signature
    assert(dark.warmth <= bright.warmth);
    assert(clean.clean >= 50); // jamais moins de debruitage que la reference

    // Reverberation : voix seche -> de-reverb hors service ; piece RT60 0.8 s (4 peignes) -> detectee, activee
    {
        auto dry = voice(0.1f, 15, 0.0003f);
        std::vector<float> wet(dry.size(), 0.f);
        for (size_t D : {1499, 1777, 2053, 2311}) { // ponytail: reverb de Schroeder minimale (peignes paralleles)
            const float g = std::pow(10.f, -3.f * D / (48000.f * 0.8f));
            std::vector<float> y(dry.size());
            for (size_t i = 0; i < y.size(); ++i) y[i] = dry[i] + (i >= D ? g * y[i - D] : 0.f);
            for (size_t i = 0; i < y.size(); ++i) wet[i] += 0.25f * y[i];
        }
        auto room = analyzeCalibration(hiss(0.0003f, 5), wet, 48000);
        std::printf("RT60 : sec %.2f s, piece %.2f s\n", clean.rt60, room.rt60);
        assert(clean.dereverbRt60 == 0 && clean.reverb == Light::Green);
        assert(room.dereverbRt60 > 0 && room.reverb != Light::Green);
        assert(room.rt60 > 0.3f && room.rt60 < 1.2f && clean.rt60 < 0.25f); // sous-evalue, mais bien separe
    }

    // Capture : downmix stereo, s'arrete a la taille demandee
    {
        Capture c; c.prepare(1000);
        std::vector<float> l(600, 1.f), r(600, 0.f);
        const float* ch[2] = {l.data(), r.data()};
        c.push(ch, 2, 600); // inactif : ignore
        c.start(800);
        c.push(ch, 2, 600);
        assert(!c.done());
        c.push(ch, 2, 600);
        assert(c.done() && c.take().size() == 800 && c.take()[0] == 0.5f);
    }

    std::puts("test_calibration OK");
    return 0;
}

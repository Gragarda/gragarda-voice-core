// Auto-test du pipeline hors temps reel (web / mobile) : plops (spec 5.5), loudness, true peak.
// test_offline [--wav] : --wav ecrit les fichiers de test et leur version traitee dans le dossier courant.
#undef NDEBUG // asserts actifs meme en Release
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include "gv/Dsp.h"
#include "gv/Pipeline.h"
#include "gv/PopRepair.h"
#include "gv/Wav.h"

using namespace gv;
static const double kTau = 6.283185307179586;
static bool gWav = false;

// "Voix" synthetique : harmoniques, f0 qui varie, syllabes a 3 Hz, souffle leger. Source continue :
// voice(seed)[i] est defini pour tout i, on peut donc y "decouper" des plans.
static float voiceAt(size_t i, float amp = 0.2f)
{
    // phase integree : f0 = 140 + 30 sin(2 pi 0.7 t) Hz, sans derive de frequence avec t
    const double t = i / 48000.0, ph = kTau * (140 * t - 30 / (kTau * 0.7) * std::cos(kTau * 0.7 * t));
    const double env = std::fmax(0, std::sin(kTau * 3 * t));
    double v = 0;
    for (int h = 1; h < 20; ++h) v += std::sin(ph * h + h) / h;
    unsigned s = (unsigned)i * 2654435761u;
    s ^= s >> 13; s *= 1274126177u;
    return (float)(amp * env * v) + 0.002f * ((s >> 8) / 8388608.f - 1.f);
}

static AudioBuffer mono(std::vector<float> x)
{
    AudioBuffer b;
    b.ch.push_back(std::move(x));
    return b;
}

// Montage sans fondu : apres la coupe k, la source reprend ~0,34 s plus loin (meme syllabe, autre phase),
// decalage ajuste pour que la forme d'onde saute d'au moins 0,1 (un vrai plop, pas un raccord chanceux)
static AudioBuffer splice(const std::vector<double>& cuts, double len, float dc)
{
    std::vector<size_t> off{0};
    for (double t : cuts) {
        const size_t c = (size_t)(t * 48000);
        size_t o = off.back() + 16373;
        while (std::fabs(voiceAt(c + o) - voiceAt(c + off.back())) < 0.1f) ++o;
        off.push_back(o);
    }
    std::vector<float> x((size_t)(len * 48000));
    size_t k = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        while (k < cuts.size() && i >= (size_t)(cuts[k] * 48000)) ++k;
        x[i] = voiceAt(i + off[k]) + (k % 2 ? dc : 0.0f);
    }
    return mono(x);
}

// Pire rapport (saut a +-2 ms de la coupe) / (plus grand saut naturel entre 2 et 12 ms autour)
static float maxResidual(const AudioBuffer& b, const std::vector<double>& cuts)
{
    float worst = 0;
    for (double t : cuts) {
        const size_t c = (size_t)(t * 48000);
        float at = 0, around = 1e-4f;
        for (size_t i = c - 576; i < c + 576; ++i)
            (i + 96 >= c && i <= c + 96 ? at : around) = std::fmax(i + 96 >= c && i <= c + 96 ? at : around, jumpAt(b, i));
        worst = std::fmax(worst, at / around);
    }
    return worst;
}

static void save(const char* name, const AudioBuffer& b)
{
    if (gWav) writeWav(name, b);
}

static std::string readFile(const char* p)
{
    std::ifstream f(p);
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}

int main(int argc, char** argv)
{
    std::setvbuf(stdout, nullptr, _IONBF, 0);
    std::string presetDir = "presets";
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--wav")) gWav = true;
        else presetDir = argv[i];
    }
    // Coupes au sommet des syllabes (sin(6 pi t) = 1) : la voix y est forte, le saut maximal
    std::vector<double> cuts;
    for (int k = 0; k < 5; ++k) cuts.push_back(1.0 / 12 + (2 * k + 1) / 3.0);

    // 1. Voix coupee/recollee a 5 endroits : 5 plops, en mode scan seul comme en mode guide
    for (int guided = 0; guided < 2; ++guided) {
        auto b = splice(cuts, 4.5, 0.0f);
        if (!guided) save("t1_splices.wav", b);
        auto pops = detectPops(b, guided ? cuts : std::vector<double>{});
        std::printf("splices %s : %zu plops\n", guided ? "guide" : "scan", pops.size());
        for (const Pop& p : pops) std::printf("  %.4f s %s x%.1f\n",(double)p.pos / 48000, toString(p.type), p.severity);
        assert(pops.size() == 5);
        for (size_t k = 0; k < 5; ++k)
            assert(std::fabs((double)pops[k].pos / 48000 - cuts[k]) < 0.002);
        const float before = maxResidual(b, cuts);
        repairPops(b, pops);
        const float after = maxResidual(b, cuts);
        std::printf("  saut/ref avant %.1f apres %.1f\n", before, after);
        assert(after < 1.5f && before > 5.0f);
        assert(detectPops(b, {}).empty());
        if (!guided) save("t1_splices_repaired.wav", b);
    }

    // 2. Meme chose avec un saut d'offset DC volontaire (+-0,05)
    {
        auto b = splice(cuts, 4.5, 0.05f);
        save("t2_dc.wav", b);
        auto pops = detectPops(b, cuts);
        std::printf("dc : %zu plops\n", pops.size());
        assert(pops.size() == 5);
        for (const Pop& p : pops) assert(p.type == PopType::Both);
        repairPops(b, pops);
        std::printf("  saut/ref apres %.1f\n", maxResidual(b, cuts));
        assert(maxResidual(b, cuts) < 1.5f);
        save("t2_dc_repaired.wav", b);
    }

    // 3. Coupe en plein silence : souffle a -60 dBFS puis saut de niveau 0,03
    {
        std::vector<float> x(96000);
        for (size_t i = 0; i < x.size(); ++i) x[i] = voiceAt(i) * 0.0f + 0.001f * std::sin(i * 0.37f) + (i >= 48000 ? 0.03f : 0.0f);
        auto b = mono(x);
        auto pops = detectPops(b, {1.0});
        std::printf("silence : %zu plops (%s)\n", pops.size(), pops.empty() ? "-" : toString(pops[0].type));
        assert(pops.size() == 1 && pops[0].type == PopType::Silence);
        repairPops(b, pops);
        assert(maxResidual(b, {1.0}) < 1.5f);
    }

    // 4. Phrase riche en "p" et "t", sans aucune coupe : ZERO plop (anti faux positif)
    {
        std::vector<float> x(48000 * 6);
        for (size_t i = 0; i < x.size(); ++i) x[i] = voiceAt(i);
        unsigned s = 7;
        double lp = 0;
        for (int k = 0; k < 18; ++k) { // une explosive au debut de chaque syllabe (attaque 0,5 ms)
            const size_t on = (size_t)((k / 3.0 + 0.002) * 48000);
            const bool p = k % 2 == 0; // "p" sourd (grave) / "t" claque (aigu)
            for (size_t i = 0; i < 4800 && on + i < x.size(); ++i) {
                s = s * 1664525u + 1013904223u;
                const double nz = (s >> 8) / 8388608.0 - 1.0;
                lp += (p ? 0.15 : 0.7) * (nz - lp);
                const double env = std::fmin(1.0, i / 24.0) * std::exp(-(double)i / (p ? 480 : 960)); // "p" ~10 ms, "t" aspire ~20 ms;
                x[on + i] += (float)(0.35 * env * (p ? lp * 3 : nz - lp));
            }
        }
        auto b = mono(x);
        save("t4_plosives.wav", b);
        auto pops = detectPops(b, {});
        std::printf("explosives sans coupe : %zu plops\n", pops.size());
        assert(pops.empty());
    }

    // 5. Pipeline complet sur chaque preset : -14 LUFS (+-0,5), <= -1 dBTP, plops comptes
    for (const char* name : {"face-camera", "voiture", "dehors", "a-plusieurs", "piece-qui-resonne"}) {
        const std::string file = readFile((presetDir + "/" + name + ".json").c_str());
        assert(!file.empty());
        // RNNoise efface en partie la voix synthetique (ce n'est pas de la parole) : on teste la loudness sans lui.
        // La 1re occurrence d'une cle gagne (parsePreset), d'ou l'ajout en tete.
        const std::string json = "{\"denoiseMin\": 0, " + file.substr(file.find('{') + 1);
        auto b = splice(cuts, 4.5, 0.0f);
        for (float& v : b.ch[0]) v *= 0.3f; // voix trop basse
        b.ch.push_back(b.ch[0]);              // stereo
        const Report r = runPipeline(b, cuts, json);
        std::printf("%s : %.2f LUFS, %.2f dBTP, %zu plops\n", name, r.process.outLufs, r.process.outTruePeakDb,
                    r.pops.size());
        assert(r.preset == name);
        assert(r.pops.size() == 5);
        assert(std::fabs(r.process.outLufs + 14.0f) < 0.5f);
        assert(r.process.outTruePeakDb <= -1.0f);
        if (!std::strcmp(name, "face-camera")) { std::printf("%s\n", r.toJson().c_str()); save("t5_processed.wav", b); }
    }

    // 6. Limiteur seul sur un signal tres chaud : true peak mesure <= plafond
    {
        std::vector<float> x(48000);
        for (size_t i = 0; i < x.size(); ++i) x[i] = 1.5f * (float)std::sin(kTau * 11025.3 * i / 48000);
        auto b = mono(x);
        limitTruePeak(b, dbToLin(-1.1f));
        std::printf("limiteur : %.2f dBTP\n", linToDb(truePeak(b)));
        assert(linToDb(truePeak(b)) <= -1.0f);
    }
    std::printf("OK\n");
    return 0;
}

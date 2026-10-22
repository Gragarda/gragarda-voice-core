// Auto-test du coeur DSP, sans OBS. Build : voir README.
#undef NDEBUG // asserts actifs meme en Release
#include "VoiceChain.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

static std::vector<float> sine(double hz, float amp, size_t n)
{
    std::vector<float> v(n);
    for (size_t i = 0; i < n; ++i)
        v[i] = amp * (float)std::sin(2 * 3.14159265358979 * hz * i / 48000.0);
    return v;
}

static float peakOf(const std::vector<float>& v, size_t from)
{
    float p = 0;
    for (size_t i = from; i < v.size(); ++i) p = std::fmax(p, std::fabs(v[i]));
    return p;
}

// Signal "voix" synthetique : harmoniques avec f0 qui varie, syllabes a 3 Hz
static std::vector<float> voice(float amp, size_t n, unsigned seed = 3)
{
    std::vector<float> x(n);
    for (size_t i = 0; i < n; ++i) {
        double t = i / 48000.0, f0 = 140 + 30 * std::sin(2 * 3.14159265 * 0.7 * t);
        double env = std::fmax(0, std::sin(2 * 3.14159265 * 3 * t)), vv = 0;
        for (int h = 1; h < 20; ++h) vv += std::sin(2 * 3.14159265 * f0 * h * t + h) / h;
        seed = seed * 1664525u + 1013904223u;
        x[i] = (float)(amp * env * vv) + 0.002f * ((seed >> 8) / 8388608.f - 1.f);
    }
    return x;
}

static double rmsDb(const std::vector<float>& v, size_t from, size_t to)
{
    double e = 0;
    for (size_t i = from; i < to; ++i) e += v[i] * v[i];
    return 10 * std::log10(e / (to - from) + 1e-20);
}

static const char* gModel = nullptr; // chemin du modele DeepFilterNet (argv[1])

int main(int argc, char** argv)
{
    const size_t N = 48000;
    if (argc > 1) gModel = argv[1];
    std::setvbuf(stdout, nullptr, _IONBF, 0); // voir les mesures meme si un assert tombe

    // Latence : debruitage 2x480 + limiteur 72 = 21.5 ms, identique debruitage on/off ;
    // + 480 (de-reverb) seulement quand il tourne
    VoiceChain vc;
    vc.prepare(48000, 2);
    assert(vc.latencySamples() == 1032);
    {
        VoiceChain v; v.prepare(48000, 1);
        VoiceChain::Macros m; m.rt60 = 0.8f; v.setMacros(m);
        assert(v.latencySamples() == 1512);
        VoiceChain::Stages st; st.dereverb = false; v.setStages(st);
        assert(v.latencySamples() == 1032);
        // l'impulsion ressort bien a la latence annoncee, de-reverb coupe puis remis
        for (bool on : {false, true}) {
            st.dereverb = on; v.setStages(st);
            std::vector<float> x(6000, 0.f); x[100] = 0.5f;
            float* p[1] = {x.data()}; v.process(p, 1, x.size());
            size_t k = 0;
            for (size_t i = 0; i < x.size(); ++i) if (std::fabs(x[i]) > std::fabs(x[k])) k = i;
            assert(k >= 100 + (size_t)v.latencySamples() && k <= 103 + (size_t)v.latencySamples());
        }
    }

    // Silence reste silence
    {
        std::vector<float> l(N, 0.f), r(N, 0.f);
        float* ch[2] = {l.data(), r.data()};
        vc.process(ch, 2, N);
        assert(peakOf(l, 0) == 0.f && peakOf(r, 0) == 0.f);
    }

    // Signal tres fort (+12 dB de trim sur un sinus a 0.9) : jamais au-dessus du plafond -1 dBFS
    {
        VoiceChain v; v.prepare(48000, 2); v.setGainDb(12);
        auto l = sine(1000, 0.9f, N), r = sine(3000, 0.9f, N);
        float* ch[2] = {l.data(), r.data()};
        v.process(ch, 2, N);
        const float ceil = std::pow(10.f, -1.f / 20.f);
        assert(peakOf(l, 0) <= ceil + 1e-5f && peakOf(r, 0) <= ceil + 1e-5f);
    }

    // Taille de bloc variable = meme sortie qu'un seul gros bloc
    {
        VoiceChain a, b; a.prepare(48000, 1); b.prepare(48000, 1);
        auto x = sine(440, 0.5f, N), y = x;
        float* pa[1] = {x.data()};
        a.process(pa, 1, N);
        size_t off = 0, sizes[] = {1, 7, 480, 1024, 333};
        for (int k = 0; off < N; k = (k + 1) % 5) {
            size_t n = std::min(sizes[k], N - off);
            float* pb[1] = {y.data() + off};
            b.process(pb, 1, n);
            off += n;
        }
        for (size_t i = 0; i < N; ++i) assert(x[i] == y[i]);
    }

    // Passe-haut : 30 Hz fortement attenue, 1 kHz passe
    {
        VoiceChain v; v.prepare(48000, 1);
        auto lo = sine(30, 0.05f, N), mid = sine(1000, 0.05f, N);
        float* p1[1] = {lo.data()}; v.process(p1, 1, N);
        VoiceChain w; w.prepare(48000, 1);
        float* p2[1] = {mid.data()}; w.process(p2, 1, N);
        assert(peakOf(lo, N / 2) < 0.3f * peakOf(mid, N / 2));
    }

    // Avant/Apres : brut retarde de la latence exacte, au volume du traite
    {
        const size_t M = 48000 * 20; // leveler stabilise (il monte a 2 dB/s)
        VoiceChain wet, ab; wet.prepare(48000, 1); ab.prepare(48000, 1); ab.setBypass(true);
        std::vector<float> x(M);
        for (size_t i = 0; i < M; ++i) // 300 Hz en syllabes, faible (-35 dB) : le traite sera bien plus fort
            x[i] = 0.03f * (float)(std::fmax(0, std::sin(2 * 3.14159265 * 3 * i / 48000.0)) * std::sin(2 * 3.14159265 * 300 * i / 48000.0));
        auto a = x, b = x;
        float* pa[1] = {a.data()}; wet.process(pa, 1, M);
        float* pb[1] = {b.data()}; ab.process(pb, 1, M);
        const double dWet = rmsDb(a, M - 48000 * 2, M), dAb = rmsDb(b, M - 48000 * 2, M);
        std::printf("Avant/Apres : traite %.1f dB, brut egalise %.1f dB\n", dWet, dAb);
        assert(std::fabs(dWet - dAb) < 1.5);
        // alignement : correlation max a 0 entre brut egalise et entree retardee
        const int L = ab.latencySamples();
        double c0 = 0, c1 = 0;
        for (size_t i = M / 2; i < M; ++i) { c0 += b[i] * x[i - L]; c1 += b[i] * x[i - L + 40]; }
        assert(c0 > c1);
    }

    // Pas de NaN sur bruit plein-echelle
    {
        VoiceChain v; v.prepare(48000, 2);
        std::vector<float> l(N), r(N);
        unsigned s = 1;
        for (size_t i = 0; i < N; ++i) { s = s * 1664525u + 1013904223u; l[i] = r[i] = (s >> 8) / 8388608.f - 1.f; }
        float* ch[2] = {l.data(), r.data()}; v.process(ch, 2, N);
        for (size_t i = 0; i < N; ++i) assert(std::isfinite(l[i]));
    }

    // Debruitage : souffle de fond (~-55 dBFS) fortement reduit. ponytail: bruit blanc synthetique, RNNoise y est erratique a fort niveau ; le vrai juge = l'ecoute
    {
        std::vector<float> a(N * 3), b;
        unsigned s = 7;
        for (auto& v : a) { s = s * 1664525u + 1013904223u; v = 0.002f * ((s >> 8) / 8388608.f - 1.f); }
        b = a;
        VoiceChain on, off; on.prepare(48000, 1); off.prepare(48000, 1); off.setDenoise(false);
        float* pa[1] = {a.data()}; on.process(pa, 1, a.size());
        float* pb[1] = {b.data()}; off.process(pb, 1, b.size());
        double ea = 0, eb = 0;
        for (size_t i = N * 2; i < N * 3; ++i) { ea += a[i] * a[i]; eb += b[i] * b[i]; }
        std::printf("souffle : %.1f dB\n", 10 * std::log10(ea / eb));
        assert(ea < eb * 0.25); // au moins -6 dB
    }

    // Debruitage off : retard identique (impulsion ressort a latencySamples)
    {
        VoiceChain v; v.prepare(48000, 1); v.setDenoise(false);
        std::vector<float> x(4000, 0.f); x[100] = 0.5f;
        float* p[1] = {x.data()}; v.process(p, 1, x.size());
        size_t k = 0;
        for (size_t i = 0; i < x.size(); ++i) if (std::fabs(x[i]) > std::fabs(x[k])) k = i;
        assert(k >= 100 + (size_t)v.latencySamples() && k <= 103 + (size_t)v.latencySamples()); // filtres = quelques echantillons d'etalement
    }

    // Qualite max (DeepFilterNet) : latence 4x480 + 72 = 41.5 ms (sans de-reverb) ; vraie voix-like ressort exactement a cette latence
    if (gModel) {
        VoiceChain v; v.prepare(48000, 1, gModel); v.setEngine(VoiceChain::Max);
        assert(v.maxAvailable());
        assert(v.latencySamples() == 1992);
        const size_t M = 48000 * 3;
        std::vector<float> x(M), y;
        unsigned s = 3;
        for (size_t i = 0; i < M; ++i) {
            double t = i / 48000.0, f0 = 140 + 30 * std::sin(2 * 3.14159265 * 0.7 * t);
            double env = std::fmax(0, std::sin(2 * 3.14159265 * 3 * t)), vv = 0;
            for (int h = 1; h < 20; ++h) vv += std::sin(2 * 3.14159265 * f0 * h * t + h) / h;
            s = s * 1664525u + 1013904223u;
            x[i] = (float)(0.05 * env * vv) + 0.002f * ((s >> 8) / 8388608.f - 1.f);
        }
        y = x;
        float* p[1] = {y.data()}; v.process(p, 1, M);
        int best = 0; double bc = -2;
        for (int d = 1900; d < 2100; ++d) {
            double c = 0, ex = 0, ey = 0;
            for (size_t i = 48000; i + 2100 < M; ++i) { c += x[i] * y[i + d]; ex += x[i] * x[i]; ey += y[i + d] * y[i + d]; }
            c /= std::sqrt(ex * ey);
            if (c > bc) { bc = c; best = d; }
        }
        std::printf("DFN : retard mesure %d (attendu ~%d), correlation %.2f\n", best, v.latencySamples(), bc);
        assert(std::abs(best - v.latencySamples()) <= 3);

        // souffle seul : DFN doit le reduire nettement
        std::vector<float> a(M), b;
        for (auto& q : a) { s = s * 1664525u + 1013904223u; q = 0.002f * ((s >> 8) / 8388608.f - 1.f); }
        b = a;
        VoiceChain on, off; on.prepare(48000, 1, gModel); off.prepare(48000, 1, gModel);
        on.setEngine(VoiceChain::Max); off.setEngine(VoiceChain::Max); off.setDenoise(false);
        float* pa[1] = {a.data()}; on.process(pa, 1, M);
        float* pb[1] = {b.data()}; off.process(pb, 1, M);
        double ea = 0, eb = 0;
        for (size_t i = M / 2; i < M; ++i) { ea += a[i] * a[i]; eb += b[i] * b[i]; }
        std::printf("DFN souffle : %.1f dB\n", 10 * std::log10(ea / eb + 1e-12));
        assert(ea < eb * 0.25);
    } else {
        std::puts("(modele DeepFilterNet non fourni : tests Qualite max sautes)");
    }

    // Modele absent : pas de crash, Max retombe sur un retard pur de meme longueur
    {
        VoiceChain v; v.prepare(48000, 1, "introuvable.tar.gz"); v.setEngine(VoiceChain::Max);
        assert(!v.maxAvailable());
        std::vector<float> x(4000, 0.f); x[100] = 0.5f;
        float* p[1] = {x.data()}; v.process(p, 1, x.size());
        size_t k = 0;
        for (size_t i = 0; i < x.size(); ++i) if (std::fabs(x[i]) > std::fabs(x[k])) k = i;
        assert(k >= 100 + (size_t)v.latencySamples() && k <= 103 + (size_t)v.latencySamples());
    }

    // De-reverb : meme retard actif ou non ; queue de reverb (peigne RT60 0.8 s) nettement reduite apres les syllabes
    {
        VoiceChain::Macros m; m.rt60 = 0.8f;
        VoiceChain v; v.prepare(48000, 1); v.setMacros(m);
        std::vector<float> x(4000, 0.f); x[100] = 0.5f;
        float* p[1] = {x.data()}; v.process(p, 1, x.size());
        size_t k = 0;
        for (size_t i = 0; i < x.size(); ++i) if (std::fabs(x[i]) > std::fabs(x[k])) k = i;
        assert(k >= 100 + (size_t)v.latencySamples() && k <= 103 + (size_t)v.latencySamples());

        const size_t M = 48000 * 4, D = 2399; // syllabes de 0.15 s toutes les 0.5 s, puis reverb
        std::vector<float> dryv(M, 0.f), wet(M);
        auto src = voice(0.05f, M);
        for (size_t i = 0; i < M; ++i) if (i % 24000 < 7200) dryv[i] = src[i] * 2.f;
        const float g = std::pow(10.f, -3.f * D / (48000.f * 0.8f));
        for (size_t i = 0; i < M; ++i) wet[i] = dryv[i] + (i >= D ? g * wet[i - D] : 0.f);
        auto a = wet, b = wet;
        VoiceChain on, off; on.prepare(48000, 1); off.prepare(48000, 1);
        on.setMacros(m); on.setDenoise(false); off.setDenoise(false);
        float* pa[1] = {a.data()}; on.process(pa, 1, M);
        float* pb[1] = {b.data()}; off.process(pb, 1, M);
        // rapport queue (100..400 ms apres la syllabe) / syllabe, sur les 2 dernieres secondes
        auto ratio = [&](const std::vector<float>& y) {
            double es = 0, et = 0;
            const size_t L = (size_t)v.latencySamples();
            for (size_t i = M / 2; i + L < M; ++i) {
                const size_t ph = i % 24000;
                if (ph < 7200) es += y[i + L] * y[i + L];
                else if (ph > 12000 && ph < 21600) et += y[i + L] * y[i + L];
            }
            return 10 * std::log10(et / es);
        };
        const double rOn = ratio(a), rOff = ratio(b);
        std::printf("de-reverb : queue/syllabe %.1f dB -> %.1f dB\n", rOff, rOn);
        assert(rOn < rOff - 4);
    }

    // Options avancees : tout coupe = entree retardee de la latence, a l'identique (meme latence)
    {
        VoiceChain v; v.prepare(48000, 1); v.setDenoise(false);
        VoiceChain::Macros m; m.rt60 = 0.8f; v.setMacros(m);
        v.setStages({false, false, false, false, false, false, false});
        auto x = voice(0.5f, N), y = x;
        float* p[1] = {y.data()}; v.process(p, 1, N);
        const size_t L = (size_t)v.latencySamples();
        for (size_t i = L; i < N; ++i) assert(std::fabs(y[i] - x[i - L]) < 1e-6f);
    }

    // Macros a 50 partout = reglages par defaut (reference validee a l'ecoute)
    {
        VoiceChain a, b; a.prepare(48000, 1); b.prepare(48000, 1);
        b.setMacros(VoiceChain::Macros{});
        auto x = voice(0.05f, N), y = x;
        float* pa[1] = {x.data()}; a.process(pa, 1, N);
        float* pb[1] = {y.data()}; b.process(pb, 1, N);
        assert(x == y);
    }

    // Chaleur : plus de grave a 100 qu'a 0 ; Clarte : plus d'aigu a 100 qu'a 0
    {
        auto run = [&](VoiceChain::Macros m, double hz) {
            // debruitage off (un sinus pur n'est pas de la voix), fenetre courte : le leveler (2 dB/s) n'a pas encore bouge
            VoiceChain v; v.prepare(48000, 1); v.setMacros(m); v.setDenoise(false);
            auto x = sine(hz, 0.01f, N);
            float* p[1] = {x.data()}; v.process(p, 1, N);
            return rmsDb(x, N / 5, N / 5 + N / 10);
        };
        VoiceChain::Macros lo, hi;
        lo.warmth = 0; hi.warmth = 1;
        const double dWarm = run(hi, 120) - run(lo, 120);
        lo = hi = VoiceChain::Macros{};
        lo.clarity = 0; hi.clarity = 1;
        const double dClar = run(hi, 4000) - run(lo, 4000);
        std::printf("Chaleur 0->100 @120 Hz : %+.1f dB, Clarte 0->100 @4 kHz : %+.1f dB\n", dWarm, dClar);
        assert(dWarm > 2.5 && dClar > 3); // plateau a 120 Hz = moitie du gain a sa frequence de coupure
    }

    // Leveler : voix faible remontee ; ecarts entre destinations = ecarts de cibles.
    // RMS brut (pauses comprises, sans filtre K) ~6 dB sous la sonie LUFS : on verifie surtout les ecarts.
    double outs[3];
    for (auto d : {VoiceChain::Stream, VoiceChain::Podcast, VoiceChain::VoiceOver}) {
        VoiceChain v; v.prepare(48000, 1); v.setDenoise(false);
        VoiceChain::Macros m; m.dest = d; v.setMacros(m);
        const size_t M = 48000 * 25;
        // sinus 300 Hz en "syllabes" : facteur de crete faible, le limiteur ne fausse pas la mesure
        std::vector<float> x(M);
        for (size_t i = 0; i < M; ++i)
            x[i] = 0.15f * (float)(std::fmax(0, std::sin(2 * 3.14159265 * 3 * i / 48000.0)) * std::sin(2 * 3.14159265 * 300 * i / 48000.0));
        float* p[1] = {x.data()}; v.process(p, 1, M);
        const double out = rmsDb(x, M - 48000 * 5, M);
        std::printf("dest %d : sortie %.1f dB RMS (cible %.0f LUFS)\n", (int)d, out, VoiceChain::targetLufs(d));
        assert(std::fabs(out - VoiceChain::targetLufs(d)) < 8);
        outs[d] = out;
    }
    for (int d = 1; d < 3; ++d) {
        const double want = VoiceChain::targetLufs(VoiceChain::Stream) - VoiceChain::targetLufs((VoiceChain::Destination)d);
        assert(std::fabs((outs[0] - outs[d]) - want) < 1.0);
    }

    std::puts("test_chain OK");
    return 0;
}

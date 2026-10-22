# Gragarda Voice — DSP core

Host-independent voice processing chain used by [Gragarda Voice](https://github.com/Gragarda/gragarda-voice) (OBS filter)
and [Gragarda Voice VST](https://github.com/Gragarda/gragarda-voice-vst) (VST3 plugin).

Chain: input trim → high-pass → denoise (RNNoise or DeepFilterNet 3) → subtractive EQ → compressor → de-esser →
additive EQ → leveler (short-term LUFS target) → look-ahead limiter. Driven by 4 macros (Clean, Warmth, Clarity, Punch).
Plus `Calibration`: analyzes 5 s of silence + 15 s of speech and suggests macros / mic correction.

Constraints: planar float audio, any block size, no allocation or lock in `process()`, fixed latency
(`VoiceChain::latencySamples()`). Denoisers are trained at 48 kHz.

## Use in a project

```sh
git submodule add https://github.com/Gragarda/gragarda-voice-core core
```
```cmake
add_subdirectory(core)
target_link_libraries(my-plugin PRIVATE gv_dsp)   # ${DFN_MODEL} = DeepFilterNet model to ship
```

## Offline pipeline (web / mobile)

`gv_offline` (`offline/`, `presets/*.json`) processes a whole file instead of a live stream:
edit-click ("plop") detection and repair guided by video cut timestamps, analysis (noise, BS.1770 loudness,
true peak, dynamic range, sibilance, room decay), then high-pass → RNNoise → expander → EQ → de-esser →
compressor → -14 LUFS → true-peak limiter (-1 dBTP). Entry point: `gv::runPipeline` (`offline/include/gv/Pipeline.h`).
Used by [gragarda-voice-web](https://github.com/Gragarda/gragarda-voice-web) through WebAssembly.

`-DGV_WITH_DFN=OFF` builds without DeepFilterNet (no Rust toolchain needed; required for WebAssembly).
In that mode the Max engine is unavailable and `VoiceChain::maxAvailable()` returns false.

## Build and test alone

DeepFilterNet (libDF) must be built next to the top-level project:

```sh
git clone https://github.com/Rikorose/DeepFilterNet.git ../DeepFilterNet
cd ../DeepFilterNet && git checkout d375b2d
cargo build --release -p deep_filter --no-default-features --features capi --lib
cd -
cmake -S . -B build && cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## License

GPL-2.0-or-later. RNNoise (BSD) in `third_party/rnnoise`, DeepFilterNet (MIT/Apache-2.0): see `licenses/`.

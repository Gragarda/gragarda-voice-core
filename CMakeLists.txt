cmake_minimum_required(VERSION 3.24)
project(gragarda-voice-core VERSION 0.4.0 LANGUAGES C CXX)

# Coeur DSP de Gragarda Voice : aucune dependance a un hote (OBS, VST3...).
# Utilisation : add_subdirectory(core) puis target_link_libraries(... gv_dsp) ; ${DFN_MODEL} = modele a installer.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON) # libs statiques liees dans un module partage (Linux)

# RNNoise (BSD, copie de obs-studio/plugins/obs-filters/rnnoise, modele integre)
file(GLOB RNNOISE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/rnnoise/src/*.c)
add_library(rnnoise STATIC ${RNNOISE_SRC})
target_include_directories(rnnoise PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/rnnoise/include)
target_compile_definitions(rnnoise PRIVATE _USE_MATH_DEFINES COMPILE_OPUS)

# DeepFilterNet 3 (libDF, Rust, MIT/Apache) : clone de github.com/Rikorose/DeepFilterNet a cote du projet
# principal, puis : cargo build --release -p deep_filter --no-default-features --features capi --lib
set(DFN_DIR "${CMAKE_SOURCE_DIR}/../DeepFilterNet" CACHE PATH "Sources DeepFilterNet (libDF compilee)")
set(DFN_MODEL "${DFN_DIR}/models/DeepFilterNet3_onnx.tar.gz" CACHE FILEPATH "Modele DeepFilterNet 3")
if(WIN32)
  set(DFN_LIB_DEFAULT "${DFN_DIR}/target/release/df.lib")
  set(DFN_SYS_LIBS ws2_32 userenv ntdll bcrypt advapi32)
elseif(APPLE)
  set(DFN_LIB_DEFAULT "${DFN_DIR}/target/release/libdf.a") # universel : lipo des deux cibles (voir CI)
  set(DFN_SYS_LIBS "-framework CoreFoundation" "-framework Security" "-framework SystemConfiguration" iconv)
else()
  set(DFN_LIB_DEFAULT "${DFN_DIR}/target/release/libdf.a")
  set(DFN_SYS_LIBS dl pthread m)
endif()
set(DFN_LIB "${DFN_LIB_DEFAULT}" CACHE FILEPATH "Bibliotheque statique libDF")

add_library(gv_dsp STATIC dsp/VoiceChain.cpp dsp/Calibration.cpp)
target_include_directories(gv_dsp PUBLIC dsp)
target_link_libraries(gv_dsp PRIVATE rnnoise PUBLIC "${DFN_LIB}" ${DFN_SYS_LIBS})
if(MSVC)
  target_compile_options(gv_dsp PRIVATE /W4 /fp:fast)
else()
  target_compile_options(gv_dsp PRIVATE -Wall -Wextra)
endif()

# Auto-tests (ctest) : par defaut seulement quand le coeur est compile seul
option(GV_CORE_TESTS "Compiler les tests du coeur DSP" ${PROJECT_IS_TOP_LEVEL})
if(GV_CORE_TESTS)
  enable_testing()
  add_executable(test_chain tests/test_chain.cpp)
  target_link_libraries(test_chain PRIVATE gv_dsp)
  add_test(NAME test_chain COMMAND test_chain "${DFN_MODEL}")
  add_executable(test_calibration tests/test_calibration.cpp)
  target_link_libraries(test_calibration PRIVATE gv_dsp)
  add_test(NAME test_calibration COMMAND test_calibration)
endif()

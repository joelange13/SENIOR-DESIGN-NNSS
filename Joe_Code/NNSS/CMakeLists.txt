cmake_minimum_required(VERSION 3.20)
project(TxRxApp C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(sdr_core
  # TX (existing)
  src/tx/OpenFile.c
  src/tx/build_CRC16_table.c
  src/tx/CRC16.c
  src/tx/crc16_table.c
  src/tx/byte_to_bits.c
  src/tx/packet_to_bits.c
  src/tx/bits_to_float_bpsk_iq.c
  src/tx/rrc_generate_taps.c
  src/tx/iq_interleaved_fir_filter.c
  src/tx/tx_pulse_shape.c
  src/tx/scramble_self_sync_after_offset.c
  src/tx/manchester_encode_after_offset.c
  src/tx/debugPrints.c
  src/tx/write_cf32.c
  # RX 
  src/rx/agc_create.c
  src/rx/agc_free.c
  src/rx/agc_bootstrap.c
  src/rx/agc_apply_inplace.c
  src/rx/fll_be_create.c
  src/rx/fll_be_process.c
  src/rx/fll_be_free.c
  src/rx/pfbcs_create.c
  src/rx/pfbcs_free.c
  src/rx/pfbcs_process.c
  src/rx/costas_bpsk_create.c
  src/rx/costas_bpsk_process.c
  src/rx/costas_bpsk_free.c
  src/rx/descramble_self_sync_after_offset.c
  src/rx/manchester_decode_after_offset.c
  src/rx/flag_present.c
  src/rx/get_opt.c
)

target_include_directories(sdr_core PUBLIC include)

add_executable(tx_bpsk src/main.c)
target_link_libraries(tx_bpsk PRIVATE sdr_core)

if (MSVC)
  target_compile_options(sdr_core PRIVATE /W4)
  target_compile_options(tx_bpsk PRIVATE /W4)
else()
  target_compile_options(sdr_core PRIVATE -Wall -Wextra -Wpedantic)
  target_compile_options(tx_bpsk PRIVATE -Wall -Wextra -Wpedantic)
endif()

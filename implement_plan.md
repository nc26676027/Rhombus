# Working Plan: Encrypted 1:N Face Embedding Recognition

> Status: **confirmed v1** (2026-06-15). Owner: face-pipeline agent. Source-of-truth for sequencing and exit criteria. If a phase slips, edit this file first, code second.

---

## 0. Goal

Ship a working end-to-end demo that, given an encrypted face database `Enc(W)` and an encrypted query `Enc(q)`, returns a single-bit indicator stream of the form `m_i = 1[<w_i, q> >= t]` plus a tag-limb stream `sum_i m_i * tag_i`, **without** revealing the database embeddings, the query, or the score vector to the computing party.

### 0.1 Quantitative targets for v1

| Param | v1 target | Stretch |
|---|---|---|
| Library size K | 10^4 | 10^5 |
| Embedding dim d | 128 | 512 |
| Poly degree N | 8192 | 16384 |
| Threshold bit-width | 1 (binary indicator) | 1 |
| End-to-end online latency (single query, 4 threads) | < 30 s on M-class laptop, AVX512 off | < 10 s with HEXL on |
| Online communication | < 50 MB / query | < 20 MB |
| Correctness | score MSE < 1e-3 vs plaintext | same |

### 0.2 Non-goals for v1

- Threshold decryption / ThFHE (separate track, see Phase 6).
- Multi-user batching (Park/IRIS CCMM is a separate backend, not on critical path).
- Dynamic enrollment (DB is fixed-size and preprocessed at startup).
- L2 normalization in HE (assume inputs are pre-normalized by the embedding model).

---

## 1. Threat model and API contract (lock these before Phase 1)

### 1.1 Threat model

v1 implements **Model A (single secret key holder = user)** from the brief section 8.

```
[User]    --sk-- owns sk, generates pk, gk, rk
[DB Own]  --pk-- receives Enc(W) after preprocessing
[Compute] --pk-- receives Enc(W) and Enc(q), runs the pipeline
[User]    --sk-- receives Enc(indicator) and Enc(tag-limb) for decryption
```

**Out of scope for v1:** collusion resistance, threshold decryption, query privacy against the DB owner.

### 1.2 Functional contract

Inputs:
- `W` (K x d) plaintext int8/int16 face embeddings, L2-normalized in float, quantized before encryption.
- `q` (d,) plaintext int8/int16 query embedding, also quantized.
- `t` (float, e.g. 0.8) plaintext threshold.

Outputs:
- `Enc(indicator)` ciphertext whose `decrypt(slots) -> m_i in {0, 1}` (slot-aligned).
- `Enc(tag-limb)` ciphertext, optionally multi-limb if `tag_i` exceeds 32 bits.

Failure modes the v1 demo must NOT crash on:
- K not a multiple of `N/2`.
- d not a power of 2 (must pad).
- Empty query (early return).

### 1.3 Confirmed decisions (2026-06-15)

1. **Quantization scheme**: Global symmetric int8 for both `W` and `q`. Range `[-127, 127]`, scale = `max_abs / 127`. This sets `mat_bits = 8`, `vec_bits = 8` in HyDia engine.
2. **`tag_i` payload**: 32-bit DB integer ID. Single `multiply_plain` per indicator ciphertext. No multi-limb needed.
3. **`Enc(indicator)` format**: Slot-packed — all `m_i` in one (or G) ciphertext(s), one decrypt call per group yields the full indicator vector. Minimizes communication and decryption cost.
4. **Start phase**: Phase 0 first — verify build environment, run existing tests, create scaffolding, then proceed sequentially.

---

## 2. Current state of the repo (snapshot 2026-06-15)

Already shipped and tested:
- `scripts/build-deps.sh` pulls Microsoft SEAL @ `119dc32` (v4.1.2), which already exposes `CKKSBootstrapper`. **No backend swap required.**
- `src/seal_api.{h,cpp}` (about 2.4k lines): `PackRLWEs`, `Expand`, custom `EcdRole`/`DcdRole`, `AuxParms`, threaded primitives.
- `src/matvec.{h,cpp}` / `src/matmul.{h,cpp}`: `RhombusMatVec`, `RhombusMatMul`. **Caveat: every public entry is `...ToSS` and uses `multiply_plain_ntt`. The matrix side is plaintext.** See `src/matvec.h:230,356,367,424` and `src/matmul.cpp:91`.
- `src/CtSerialization.{h,cpp}`: byte/string-level serialization of `Ciphertext`, `PublicKey`, `GaloisKeys`. No socket code.
- `src/matrix.h`: `MatMeta`, `LargeMatMeta`, `SPP_PackRLWEs`, `SPP_Expand`.
- `src/ring_op.h`: int8/16/32/64 random vectors, real vectors.
- `test/test_matvec.cpp` and `test/test_matmul.cpp` exercise PP-MVM only; no slot-domain scoring test.
- Default params are `N=8192, mod_bits=37, coeff_mod_bits={50,50,60}` (3 CKKS levels).

Not present in repo:
- Any path that takes two ciphertexts and produces a slot-aligned output.
- Threshold polynomial, indicator cleaning, indicator*tag.
- Bootstrap, CoeffToSlot, ThFHE/ThresholdDecrypt.
- Network transport.

---

## 3. Gaps mapped to the brief

| Brief section | Status | Phase that closes it |
|---|---|---|
| Section 6.2 HyDia-style diagonal MVM | missing | P1 |
| Section 3 Rhombus-coeff EE-MVM baseline | missing (would need patch, not required) | P1.5 (deferred) |
| Section 6.2 threshold polynomial | missing | P2 |
| Section 6.2 cleaning / BLEACH | missing | P2 |
| Section 6.2 indicator*tag | missing | P3 |
| Section 6.2 count / total | missing | P3 |
| Network (brief 2.3 says excluded) | missing | P5 |
| Threshold decryption (section 8 model C) | missing | P6 (separate) |
| Bootstrap | missing | P4 (conditional) |

---

## 4. Phased plan

Each phase has a single owner exit criterion. Do not start phase N+1 until N's exit criterion passes.

### Phase 0 — Preparation (0.5 day)

Deliverable: buildable repo + smoke test.

- [ ] Verify `bash scripts/build-deps.sh && bash scripts/build.sh` produces `build/bin/matvec` and `build/bin/matmul` on a Linux x86-64 box. Apple silicon is allowed to skip HEXL; do not block.
- [ ] Re-run the two existing test binaries; capture baseline numbers.
- [ ] Add a `scripts/run-tests.sh` that wraps the two binaries and exits non-zero on any `max_diff != 0`.
- [ ] Pin v1 quantization in `src/face/quant.h`: global symmetric int8, `QUANT_BITS = 8`, `QUANT_RANGE = 127`, `mat_bits = 8`, `vec_bits = 8`. *(Decision confirmed 2026-06-15)*
- [ ] Create `src/face/` directory structure and `src/face/CMakeLists.txt`.
- [ ] Add `RHOMBUS_BUILD_FACE` option to top-level `CMakeLists.txt` (default OFF).
- [ ] Decide whether to enable `USE_HEXL=ON` and `USE_ZSTD=ON` for the dev box; if yes, re-run baselines with them on.

Exit: `scripts/run-tests.sh` is green; baseline numbers recorded in this file.

### Phase 1 — HyDia engine (slot-domain EE-MVM) (2–3 days)

Goal: `Enc(W) * Enc(q) -> Enc(score)` where `score_i` lives in CKKS **slots**, ready for threshold.

New files:
- `src/face/hydia.h` / `src/face/hydia.cpp` — class `HyDiaEngine`
- `src/face/diagonal_layout.h` / `.cpp` — diagonal-block packing helpers
- `test/test_hydia.cpp` — correctness + per-group latency

API surface (rough):

```cpp
class HyDiaEngine {
public:
  HyDiaEngine(size_t d, size_t K, uint32_t mat_bits, uint32_t vec_bits, uint32_t t_bits);
  void load_db(const int8_t* W /* K*d */);
  vector<Ciphertext> encrypt_db();             // Enc(DB[g,h]) for all g,h
  Ciphertext encrypt_query(const int8_t* q);    // Enc(q | q | ... | q) repeated packing
  vector<Ciphertext> matvec(const vector<Ciphertext>& enc_db, const Ciphertext& enc_q);
private:
  // delayed relin/rescale pattern from brief section 3.4
};
```

Algorithm (from brief 2.1 + 6.2):
1. Split `W` into `G = ceil(K / (N/2))` groups of `d` rows.
2. For each group, store `d` "diagonal" ciphertexts: `DB[g, h] = Enc(Diag_h(W_g))`, each with `N/2` packed diagonals.
3. Pre-rotate query: cache `Rot(enc_q, h)` for `h = 0 .. d-1`. Each `Rot` is one `apply_galois`.
4. Per group online: `acc_deg2 = sum_h DB[g,h] * Rot(enc_q, h)`, **single** relin + rescale per group. This is the key "delayed relin" pattern that keeps the level count to 1.
5. Output: `G` ciphertexts, each holding `N/2` scores.

Implementation notes:
- Use the existing `seal::Evaluator::multiply` from SEAL directly; **do not** route through `multiply_plain_ntt`.
- Use a per-group `sum` accumulator (brief 3.4) to avoid the `i1_count * tau_count` relin explosion.
- Set `remain_n_mod = 2` initially; drop to `1` only if a noise-budget profile says so.
- `galois_keys` for `Rot` must be generated for steps 1, 2, 4, ..., `N/2`. `CtSerialization::GenGaloisKey` already does this if you pass a custom `galois_elts` list (see `src/CtSerialization.cpp:48-58`).

Test:
- K=2048, d=128, N=8192, `max_diff < 1` against plaintext reference at 8-bit quantization.
- Latency: per-group online `multiply + rotate` chain under 50 ms (single thread, no HEXL).

Exit: `build/bin/test_hydia` shows 0 max diff on K=2048 d=128 and prints G per-group online latency.

### Phase 2 — Threshold + indicator cleaning (1–2 days)

Goal: given `Enc(score)`, output `Enc(indicator)` with `indicator_i` close to `{0, 1}`.

New files:
- `src/face/threshold.h` / `.cpp` — `ThresholdPoly` (Chebyshev + minimax)
- `src/face/clean.h` / `.cpp` — `CleanIndicator` (squaring + square-root, brief 6.2 step "cleaning")

Chebyshev/minimax fit:
- Domain: `[0, 1]` (or scaled by 1.0 for cos sim, since `<w,q> = cos(theta)` after L2 norm).
- Threshold: `t = 0.8`.
- Target: `1[x >= t]`, approximated as degree `D` polynomial `P_D(x - t)`.
- Constraint: `D` must fit the remaining CKKS level budget. With 3 levels in default params and 1 used by HyDia, we have 2 levels left = 4 mults available. So pick `D <= 7` and a `clean_factor` of 2.
- Place coefficients in a `Plaintext` (NTT form, scale = current ciphertext scale).

Cleaning (BLEACH-style):
- Apply `P_D` once -> `~m_i in [0,1]`.
- Square twice -> concentrates near 0 and 1, then threshold-square-root trick. Or just use `clean_factor = 2` repeated squaring.
- Each squaring is one `multiply` + `relinearize` + `rescale_to_next` -> eats 1 level. With 2 levels left after P_D, max 2 squarings.

API:

```cpp
class ThresholdPoly {
public:
  ThresholdPoly(double t, uint32_t degree, double scale);
  Plaintext coeff_plaintext(const SEALContext& ctx) const;
};
Plaintext clean_coefficients(uint32_t rounds);  // returns the "amplification" plaintext
```

Test:
- Inputs: `Enc(score)`, scale 2^40. Output: `decrypt(Enc(indicator))` should be within `1e-2` of `{0,1}`.
- Check that `sum_i m_i` (hom count) equals the true count to within plus or minus 1.

Exit: `build/bin/test_threshold` shows clean `{0,1}` indicators on a synthetic K=2048 score set.

### Phase 3 — Indicator * tag + count (0.5 day)

Goal: output `Enc(tag-limb)` and `Enc(count)`.

New files:
- `src/face/tagselect.h` / `.cpp` — `IndicatorTagSelect` (plaintext tag vector * ciphertext indicator)
- `src/face/count.h` / `.cpp` — `HomCount` (sum of slot values)

Implementation:
- `tagselect`: pack `tag_i` (32-bit int DB ID, confirmed) as a plaintext (one per slot, in NTT form), `multiply_plain` against the indicator ciphertext. Single level consumed. No multi-limb needed.
- `count`: rotate-and-sum to fold all slots into slot 0, then `rescale`. Or just do `multiply_plain` with a `[1, 1, ..., 1]` mask followed by `rotate-and-add`. Each `rotate` is one Galois.

Test:
- Run on K=1024 with random 32-bit tags. Decrypt and compare tag-limb value with the true sum.

Exit: `build/bin/test_tagselect` returns correct `count` and `tag` on random data.

### Phase 4 — Bootstrap integration (conditional, 1 day)

Trigger: Phase 1 + 2 + 3 together consume more than 3 CKKS levels (the default). At d=128 with `D=7` and `clean_factor=2`, profile should give the answer. If levels are tight:

- [ ] Bump default `coeff_mod_bits` to `{60, 60, 60, 60}` (4 levels + 1 special modulus for bootstrap).
- [ ] In `HyDiaEngine`, expose `set_bootstrap_point(LEVEL)` so caller can decide where to refresh.
- [ ] Add `SEAL CKKSBootstrapper` instance, generate `BSKEY` (SEAL `GaloisKeys` style keys for bootstrap), serialize through `CtSerialization`.
- [ ] Insert `bootstrapper.bootstrap()` after `ThresholdPoly` evaluation but before cleaning if needed.

NOT to do: do not change the SEAL fork, do not bump past v4.1.2 unless a specific bug blocks.

Exit: end-to-end pipeline runs at d=128 with 0 max diff and a recorded latency number.

### Phase 5 — Network transport (1 day)

Goal: replace the in-process `string`-based handshakes in `test_matvec.cpp:53-72` with a real socket pair.

New files:
- `src/face/net/transport.h` / `.cpp` — `Sender`, `Receiver` over TCP, length-prefixed framing, no TLS in v1 (note in README).
- `src/face/net/keys_exchange.h` — wraps `CtSerialization` with a `send_pk / send_gk / send_rk / send_ct` interface.

Decisions to make:
- Library: zeroMQ (zmtp) is a good fit; cppzmq is header-only. Asio + plain TCP is a smaller dependency footprint.
- Length-prefix: 4-byte little-endian length header before each payload.
- Reconnect: out of scope v1, fail hard on disconnect.

Test:
- Two processes on localhost: `face_server` listens, `face_client` sends `Enc(q)`, receives `Enc(indicator)` and `Enc(tag)`.
- Total bytes on wire for K=10^4, d=128 should be under 50 MB (target).

Exit: a runnable `bin/face_server` + `bin/face_client` pair that completes one round-trip.

### Phase 6 — Threshold decryption (separate, optional for v1)

Out of scope for the v1 demo. Plan only:
- Microsoft SEAL does not include ThFHE. Do not try to bolt it onto SEAL.
- Candidate libs: OpenFHE (BFV threshold; CKKS threshold in newer versions), Lattigo (Go), custom DKG over the existing BFV/CKKS code.
- Likely a 2-3 week effort; gate behind a v2 milestone and a real security review.

### Phase 7 — End-to-end benchmark + memory/storage profile (1 day)

Goal: publish numbers matching section 0.1.

- [ ] Measure online latency, communication, peak RSS for K in {10^3, 10^4}, d in {128, 256}.
- [ ] Measure DB preprocessing time (one-shot, offline).
- [ ] Compare against Rhombus PP-MVM baseline running at the same N.
- [ ] Dump results into `docs/bench-v1.md`.

Exit: `docs/bench-v1.md` exists and answers the section 0.1 table with measured numbers.

---

## 5. File/module layout proposal

```
src/
  face/                       # new namespace, gated on RHOMBUS_BUILD_FACE
    quant.h                   # W/q quantization
    hydia.{h,cpp}             # Phase 1
    diagonal_layout.{h,cpp}   # Phase 1
    threshold.{h,cpp}         # Phase 2
    clean.{h,cpp}             # Phase 2
    tagselect.{h,cpp}         # Phase 3
    count.{h,cpp}             # Phase 3
    net/
      transport.{h,cpp}       # Phase 5
      keys_exchange.{h,cpp}   # Phase 5
  face/CMakeLists.txt         # adds face/ subdir
  seal_api.{h,cpp}            # unchanged
  matvec.{h,cpp}              # unchanged (PP-MVM remains as ablation)
  matmul.{h,cpp}              # unchanged
test/
  test_hydia.cpp              # Phase 1
  test_threshold.cpp          # Phase 2
  test_tagselect.cpp          # Phase 3
  test_face_e2e.cpp           # Phase 7
scripts/
  build.sh                    # add -DRHOMBUS_BUILD_FACE=ON default
  run-tests.sh                # Phase 0
  run-face-e2e.sh             # Phase 5
docs/
  bench-v1.md                 # Phase 7
```

The `RHOMBUS_BUILD_FACE` gate keeps the existing `matvec` / `matmul` tests usable as the original-paper ablation.

---

## 6. Test strategy

- Unit tests: one binary per new module under `test/`. Each test should print `max_diff` and an explicit pass/fail. Use the existing `CompVector` style in `test_matvec.cpp:120-130` as the template.
- Determinism: every test seeds a `std::mt19937` with a fixed seed. Random vector generation already lives in `src/ring_op.h`.
- Reference comparisons: always compare against a plaintext reference run on the same `int` quantized inputs, never against SEAL's CKKS float `encode` (different rounding).
- Network test: spawn server and client as child processes from a single test harness, fail if either crashes.
- Performance: only run perf tests when `RHOMBUS_PERF_TESTS=ON`; do not gate CI on them.

---

## 7. Build & CI

- Add `RHOMBUS_BUILD_FACE` option to the top-level `CMakeLists.txt` (default `OFF`, flip to `ON` when Phase 1 starts).
- Keep `USE_HEXL=ON` and `USE_ZSTD=ON` available but do not require them for CI.
- Document the macOS skip path: HEXL is x86-64 only; `scripts/build-deps.sh` will produce a non-HEXL SEAL on Apple silicon and that's fine.
- CI: GitHub Actions `ubuntu-22.04`, `make -j`, then `build/bin/matvec` and `build/bin/matmul`. Add the new `test_hydia`, `test_threshold`, `test_tagselect` to the same runner.

---

## 8. Risks and how we handle them

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| 3 CKKS levels not enough at d=128 | medium | P1-P3 must include bootstrap | Phase 4 is gated; bump `coeff_mod_bits` rather than swap SEAL |
| `Rot(enc_q, h)` for `h = 0..d-1` blows the rotation-key size | medium | 1:N slower | Generate only the `d` rotations actually used; store via existing `CtSerialization::GenGaloisKey(custom_list)` |
| Tag multi-limb overflows level budget | ~~eliminated~~ | — | Confirmed: 32-bit int ID, single `multiply_plain`, no multi-limb |
| Diagonal block packing wastes slots for K not a power of 2 | high | bandwidth | Pad to next multiple of `N/2`; record the waste in `bench-v1.md` |
| CKKS noise from `Rot + Multiply + relin + rescale` chain at d=512 | high | wrong results | Constrain v1 to d=128; document d=512 as a stretch goal that needs Phase 4 |
| Threshold decryption requested mid-development | low | scope creep | Phase 6 is intentionally separate; refuse to merge into v1 |

---

## 9. References

- Brief: `agent_brief_encrypted_face_mvm.md` in this repo.
- Rhombus paper: https://eprint.iacr.org/2024/1611.pdf
- HyDia-style diagonal ct-ct MVM: brief 2.1 and references therein.
- Microsoft SEAL bootstrap docs: `native/src/seal/c/ckksbootstrapper.h` of `microsoft/SEAL@119dc32`.
- Existing serialization API: `src/CtSerialization.h`.

---

## 10. Definition of Done for v1

A single command (`bash scripts/run-face-e2e.sh`) starts a server, runs the client, prints:
1. `count` and `tag` decrypted at the user.
2. `max_diff` between homomorphic result and plaintext reference is `0` for `indicator` and within plus or minus 1 for `count`.
3. Online latency and communication numbers within the section 0.1 targets.
4. `docs/bench-v1.md` is updated with the run.

When all four hold, v1 is shipped and Phase 6 (ThFHE) becomes the next planning exercise.

---

## 11. Timeline estimate (added 2026-06-15)

| Phase | Estimated effort | Dependencies |
|-------|-----------------|--------------|
| Phase 0 — Preparation | 0.5 day | None |
| Phase 1 — HyDia engine | 2–3 days | Phase 0 exit |
| Phase 2 — Threshold + cleaning | 1–2 days | Phase 1 exit |
| Phase 3 — Tag select + count | 0.5 day | Phase 2 exit |
| Phase 4 — Bootstrap (conditional) | 1 day | Level-budget profile from P1+P2+P3 |
| Phase 5 — Network transport | 1 day | Phase 3 exit |
| Phase 7 — Benchmark | 1 day | Phase 5 exit |
| **Total (core P0–P3)** | **4–7 days** | — |
| **Total (full v1 P0–P7)** | **7–10 days** | — |

Development starts from Phase 0, strictly sequential. Phase 4 is triggered only if level-budget profiling shows it is needed.

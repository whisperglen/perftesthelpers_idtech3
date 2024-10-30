I used this harness to test my work for quake3, when converting altivec optimisations to SSE and neon.

Binary dumps can be found in the sister repo: https://github.com/whisperglen/perftesthelpers_testdata.

Commited so far:<br>
-rsqrt conversion
-sndmix conversion (S_WriteLinearBlastStereo16)

Results:<br>
| function          |  x86   |  x64   |
| :---------------- | :----: | :----: |
| rsqrt_math        |  1     | 1      |
| rsqrt_q3          |  0.64  | 0.98   |
| rsqrt_sse_approx  |  0.54  | 0.80   |
| rsqrt_sse         |  0.40  | 0.42   |
| sndmix_scalar     |  1     | 1      |
| rsqrt_mmx         |  0.22  | n/a    |
| rsqrt_sseasm      |  0.13  | n/a    |
| rsqrt_sse         |  0.14  | 0.16   |

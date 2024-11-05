I used this harness to test my work for quake3, when converting altivec optimisations to SSE and neon.

Binary dumps can be found in the sister repo: https://github.com/whisperglen/perftesthelpers_testdata.

Commited so far:<br>
-rsqrt conversion<br>
-sndmix conversion (S_WriteLinearBlastStereo16 and S_PaintChannelFrom16)<br>
-dotproduct experiments (SSE4.1 is good enough, SSE3 improves when run in a row on the same dataset)<br>
-diffusecolor<br>
-lerpmeshvertexes<br>
-projectdlighttexture<br>

These results are from a destop with an i7-7700. When running on a laptop i7-6500u the SSE code performs a little better.

Results:<br>
| function          |  x86   |  x64   | x64/x86|
| :---------------- | :----: | :----: | :----: |
| rsqrt_math        |  1     | 1      | 0.67   |
| rsqrt_q3          |  0.64  | 0.98   | 0.65   |
| rsqrt_sse_precise |  0.54  | 0.80   | 0.53   |
| rsqrt_sse         |  0.40  | 0.42   | 0.28   |
|                   |        |        |        |
| sndmix_scalar     |  1     | 1      | 0.67   |
| sndmix_mmx        |  0.22  | n/a    | n/a    |
| sndmix_sseasm     |  0.13  | n/a    | n/a    |
| sndmix_sse        |  0.14  | 0.16   | 0.10   |
|                   |        |        |        |
| dotprod           |  1     | 1      | 1.1    |
| dotprod_sse       |  1.1   | 0.96   | 1.05   |
| dotprod_sse_dp    |  0.64  | 0.86   | 1.03   |
|                   |        |        |        |
| diffusecolor      |  1     | 1      | 0.95   |
| diffusecolor_sse  |  0.4   | 0.38   | 0.36   |
|                   |        |        |        |
| lerpmesh          |  1     | 1      | 0.82   |
| lerpmesh_sse      |  0.82  | 0.79   | 0.62   |
|                   |        |        |        |
| projectdlight     |  1     | 1      | 0.94   |
| projectdlight_sse |  1.02  | 1.02   | 0.96   |
|                   |        |        |        |
| sndpaint          |  1     | 1      | 0.91   |
| sndpaint_sse      |  0.31  | 0.32   | 0.29   |

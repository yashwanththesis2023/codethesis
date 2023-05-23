\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  y
Subject To
 c1: x1 + 27 r1 <= 15
 c2: - y + y1 <= -28
 c3: x2 + 19 r2 <= 15
 c4: - y + y2 <= -20
 c5: x3 + 2.147483648e+09 r3 <= -2.147483632e+09
 c6: - y + y3 <= 0
 c7: x1 + 27 r1 - 16 p12 - 16 q12 - x2 <= -1
 c8: - 27 r1 + y1 - 16 p12 + 16 q12 - y2 <= -12
 c9: x1 - 16 p12 + 16 q12 - x2 - 19 r2 >= -15
 c10: y1 - 16 p12 - 16 q12 + 19 r2 - y2 >= -12
 c11: x1 + 27 r1 - 16 p13 - 16 q13 - x3 <= -1
 c12: - 27 r1 + y1 - 16 p13 + 16 q13 - y3 <= -12
 c13: x1 - 16 p13 + 16 q13 - x3 - 2.147483648e+09 r3 >= 2.147483632e+09
 c14: y1 - 16 p13 - 16 q13 + 2.147483648e+09 r3 - y3 >= -32
 c15: x2 + 19 r2 - 16 p23 - 16 q23 - x3 <= -1
 c16: - 19 r2 + y2 - 16 p23 + 16 q23 - y3 <= -4
 c17: x2 - 16 p23 + 16 q23 - x3 - 2.147483648e+09 r3 >= 2.147483632e+09
 c18: y2 - 16 p23 - 16 q23 + 2.147483648e+09 r3 - y3 >= -32
Bounds
 r1 <= 1
 r2 <= 1
 r3 <= 1
Binaries
 p12 q12 p13 q13 p23 q23
Generals
 x1 r1 y1 x2 r2 y2 x3 r3 y3
End

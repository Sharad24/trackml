#pragma once

namespace FW {
struct Layer {
  double minr, avgr, maxr;
  double minz, avgz, maxz;
  int count;
  int type;

  double var0, var1;
};
}

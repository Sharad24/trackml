#pragma once

namespace FW {

struct triple {
  int x, y, z; //hit ids
  triple() {}
  triple(int a, int b, int c) : x(a), y(b), z(c) {}
};
bool operator<(const triple&a, const triple&b) {
  if (a.z != b.z) return a.z < b.z;
  if (a.y != b.y) return a.y < b.y;
  return a.x < b.x;
}
bool operator==(const triple&a, const triple&b) {
  return a.x==b.x && a.y==b.y && a.z==b.z;
}

}

Algorithms for high school CS

## Usage:
```cpp
int ax, ay, bx, by, cx, cy, dx, dy;
std::cin >> ax >> ay >> bx >> by >> cx >> cy >> dx >> dy;

auto i = geo::seg_intersection<int>({{ax, ay}, {bx, by}}, {{cx, cy}, {dx, dy}});

std::cout << i.value() << std::endl;
```

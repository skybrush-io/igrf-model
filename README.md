igrf-model
==========

A pure Python implementation of the IGRF-13 magnetic model, with zero
dependencies.

The implementation was ported from the original Fortran code. Being a pure
Python implementation, don't expect stellar performance, but it is good enough
for single-point queries.

Usage
-----

```python
from igrf_model import IGRFModel

model = IGRFModel.get(version=13)
vec = model.evaluate(47.498, 19.04, date=datetime(2021, 12, 22))
print(vec.declination, vec.inclination)
```

Refer to the docstrings of the `IGRFModel` and `MagneticVector` classes for
more details.

License
-------

Copyright 2021-2022 CollMot Robotics Ltd.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

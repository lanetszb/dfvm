/* MIT License
 *
 * Copyright (c) 2020 Aleksandr Zhuravlyov and Zakhar Lanets
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "Boundary.h"
#include <algorithm>
#include <iterator>

Boundary::Boundary(std::shared_ptr<Props> props,
                   std::shared_ptr<Sgrid> sgrid) :
        _props(props),
        _sgrid(sgrid) {}

void Boundary::shiftBoundaryFaces(Eigen::Ref<Eigen::VectorXui64> faces,
                                  const uint8_t &axis) {

    for (uint64_t i = 0; i < faces.size(); i++) {
        auto &face = faces[i];
        auto neighborCell = _sgrid->_neighborsCells[face][0];
        for (auto &neighborFace: _sgrid->_neighborsFaces[neighborCell]) {
            auto neighborFaceAxis = _sgrid->calculateAxisFace(neighborFace);
            if (neighborFace != face and neighborFaceAxis == axis) {
                face = neighborFace;
                break;
            }
        }
    }

}

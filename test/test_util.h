/* Copyright (C) 2017 haniu (niuhao.cn@gmail.com)
 *
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef __IEXP_TEST_UTIL__
#define __IEXP_TEST_UTIL__

////////////////////////////////////////////////////////////
// import header files
////////////////////////////////////////////////////////////

#include <common/environment.h>

#include <common/namespace.h>

#ifdef IEXP_MGL2
#include <mgl2/mgl.h>
#endif

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define __D_EQ_IN(a, b, e)                                                     \
    ((std::abs((a) - (b)) < (e)) ||                                            \
     (std::isnan((double)(a)) && std::isnan((double)(b))))

#define __D_EQ(a, b) __D_EQ_IN(a, b, std::numeric_limits<double>::epsilon())

#define __D_EQ9(a, b) __D_EQ_IN(a, b, 1e-9)

#define __D_EQ7(a, b) __D_EQ_IN(a, b, 1e-7)

#define __D_EQ6(a, b) __D_EQ_IN(a, b, 1e-6)

#define __D_EQ2(a, b) __D_EQ_IN(a, b, 1e-2)

#define __D_EQ3(a, b) __D_EQ_IN(a, b, 1e-3)

#define __D_EQ_Nep(a, b, n)                                                    \
    __D_EQ_IN(a, b, (n)*std::numeric_limits<double>::epsilon())

#define __F_EQ_IN(a, b, e)                                                     \
    ((std::abs((a) - (b)) < (e)) ||                                            \
     (std::isnan((float)(a)) && std::isnan((float)(b))))

#define __F_EQ(a, b) __F_EQ_IN(a, b, std::numeric_limits<float>::epsilon())

#define __F_EQ5(a, b) __F_EQ_IN(a, b, 0.00001)

#define __F_EQ7(a, b) __F_EQ_IN(a, b, 1e-7)

#define __F_EQ8(a, b) __F_EQ_IN(a, b, 1e-8)

// clang-format off
#define except_begin()                                                         \
    do {                                                                       \
        bool ok = false;                                                       \
        try

#define except_str(str)                                                        \
        catch (std::exception & e) {                                           \
            ok = (e.what() != nullptr) && (strcmp(e.what(), (str)) == 0);      \
            REQUIRE(ok);                                                       \
        }                                                                      \
    }                                                                          \
    while (0)
// clang-format on

////////////////////////////////////////////////////////////
// type definition
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// global variants
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
// interface declaration
////////////////////////////////////////////////////////////

IEXP_NS_END

#endif /* __IEXP_TEST_UTIL__ */

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

IEXP_NS_BEGIN

////////////////////////////////////////////////////////////
// macro definition
////////////////////////////////////////////////////////////

#define __D_EQ(a, b)                                                           \
    ((std::abs((a) - (b)) < (std::numeric_limits<double>::epsilon())) ||       \
     (std::isnan((a)) && std::isnan((b))))

#define __D_EQ_IN(a, b, e)                                                     \
    ((std::abs((a) - (b)) < (e)) || (std::isnan((a)) && std::isnan((b))))

#define __D_EQ9(a, b) __D_EQ_IN(a, b, 1e-9)

#define __D_EQ7(a, b) __D_EQ_IN(a, b, 1e-7)

#define __D_EQ6(a, b) __D_EQ_IN(a, b, 1e-6)

#define __D_EQ_Nep(a, b, n)                                                    \
    __D_EQ_IN(a, b, (n)*std::numeric_limits<double>::epsilon())

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

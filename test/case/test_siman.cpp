#include <../test/test_util.h>
#include <catch.hpp>
#include <iostream>
#include <siman/siman.h>

using namespace iexp;

namespace Eigen {

static int ctor = 0, dtor = 0;

class siman_test
{
  public:
    static void test_get_set()
    {
        siman<int> s1;

        s1.n_tries(5);
        REQUIRE(s1.m_n_tries == 5);

        s1.iters_fixed_T(7);
        REQUIRE(s1.m_iters_fixed_T == 7);

        s1.cooling(1.1, 2.2, 3.3, 4.4);
        REQUIRE(s1.m_k == 1.1);
        REQUIRE(s1.m_t_initial == 2.2);
        REQUIRE(s1.m_mu_t == 3.3);
        REQUIRE(s1.m_t_min == 4.4);
    }

    static void test_basic()
    {
        siman<int> s1;

        // energy
        s1.energy([](int x) { return x + 1; });
        std::shared_ptr<siman<int>::state> p(
            new siman<int>::state(s1, new int(1)));
        REQUIRE(s1.s_Efunc(p.get()) == 2);

        // step
        s1.step_size(3.0);
        REQUIRE(s1.m_step_size == 3.0);
        s1.step([](rand::rng &r, int &x, double step_size) {
            REQUIRE(step_size == 999);
            x = 123;
        });
        s1.s_step(nullptr, p.get(), 999);
        REQUIRE(*(int *)p->m_cfg == 123);

        // metric
        std::shared_ptr<siman<int>::state> p2(
            new siman<int>::state(s1, new int(98)));
        s1.metric([](const int &x, const int &y) -> double {
            REQUIRE(x == 123);
            REQUIRE(y == 98);
            return x - y;
        });
        double v = s1.s_metric(p.get(), p2.get());
        REQUIRE(v == (123 - 98));

        // print
        s1.print([](const int &x) { REQUIRE(x == 123); });
        s1.s_print(p.get());
    }

    static void test_copy()
    {
        siman<int> s1([](const int &s, int &d) { d = s; },
                      [](const int &s) -> int & { return *(new int(s)); });

        std::shared_ptr<siman<int>::state> p(
            new siman<int>::state(s1, new int(1)));
        std::shared_ptr<siman<int>::state> p2(
            new siman<int>::state(s1, new int(2)));
        s1.s_copy(p.get(), p2.get());
        REQUIRE(p->m_cfg != p2->m_cfg);
        REQUIRE(*(int *)p->m_cfg == 1);
        REQUIRE(*(int *)p2->m_cfg == 1);

        *(int *)p->m_cfg = 9;
        p2.reset((siman<int>::state *)s1.s_copy_construct(p.get()));
        REQUIRE(*(int *)p->m_cfg == 9);
        REQUIRE(*(int *)p2->m_cfg == 9);

        siman<int>::state *ps = new siman<int>::state(s1, new int(3));
        s1.s_destroy(ps);

        // object
        siman<to> s2([](const to &s, to &d) { d = s; },
                     [](const to &s) -> to & { return *(new to(s)); });
        std::shared_ptr<siman<to>::state> q(
            new siman<to>::state(s2, new to(1)));
        std::shared_ptr<siman<to>::state> q2(
            new siman<to>::state(s2, new to(2)));
        s2.s_copy(q.get(), q2.get());
        REQUIRE(q->m_cfg != q2->m_cfg);
        REQUIRE(*((to *)q->m_cfg)->x == 1);
        REQUIRE(*((to *)q2->m_cfg)->x == 1);
        REQUIRE(ctor == 2);

        *((to *)q->m_cfg)->x = 9;
        q2.reset((siman<to>::state *)s2.s_copy_construct(q.get()));
        REQUIRE(*((to *)q->m_cfg)->x == 9);
        REQUIRE(*((to *)q2->m_cfg)->x == 9);
        REQUIRE(ctor == 3);
    }

    struct to
    {
        to(int v)
            : x(new int(v))
        {
            ++ctor;
        }
        ~to()
        {
            ++dtor;
        }

        int *x;

        to(const to &other)
            : x(new int(*other.x))
        {
            ++ctor;
        }
        to &operator=(const to &other)
        {
            *x = *other.x;
            return *this;
        }
    };
};
}

TEST_CASE("siman")
{
    siman_test::test_get_set();
    siman_test::test_basic();

    siman_test::test_copy();
    REQUIRE(dtor == 3);
}

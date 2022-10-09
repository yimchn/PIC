#pragma once

#include <functional>
#include <iostream>

/**
 * @brief 装饰器基类，内部封装一个函数器，用于指定要被包装的函数
 *
 * @tparam RET 被包装函数的返回值
 * @tparam ARGS 被包装函数的参数
 */
template <typename RET, typename... ARGS>
class Decorator {
   protected:
    std::function<RET(ARGS...)> decorated_func;  // 被装饰的函数

   public:
    // 针对函数器的构造函数
    explicit Decorator(const std::function<RET(ARGS...)> &p_func)
        : decorated_func(p_func) {}

    // 针对函数指针的构造函数
    explicit Decorator(RET (*p_func)(ARGS...))
        : Decorator(std::function<RET(ARGS...)>(p_func)){};

    // 针对函数器重载()运算符
    virtual RET operator()(ARGS... args) = 0;
};
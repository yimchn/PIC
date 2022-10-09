#pragma once

#include <chrono>

#include "decorator.h"

/**
 * @brief 通过宏定义的匿名函数来处理函数返回值为void时编译报错的问题
 * 
 * @param content 被装饰的函数上下文（对与成员函数，可以包括类的初始化等等）
 * @param text 针对计时器的功能说明
 *
 * 当函数返回值为void时，模板替换的结果为void result =
 * this->decorated_func(args...); 由于C++中没有void类型变量，所以编译会报错。
 * 这里通过返回值为void*的lambda函数进行二次封装来解决这个问题，此处使用宏定义的形式。
 *
 * 用法：
 * TIMING(Exapmle example_class(param); example_class.example_func;,"tip");
 *
 * TIMING类详解
 * TIMING宏将Timer类初始化时的参数与返回值均设置为void*。 在Timer<void*, void*>
 * 后面有两个括号，前一个括号为Timer类的初始化，后一个括号为初始化后的函数调用，后一个括号
 * 中的内容就是被修饰函数的参数。
 * 初始化部分宏展开的结果为：
 * Timer<void* ,void*>(void*, "tip")
 * 由于content内容即为标准C++语句，可以包括类型的初始化、被装饰函数的运行等。在do后会执行，
 * 所以即便初始化时候的被装饰的函数为void*也不影响其功能的正常运行，反而因为将函数与参数写
 * 在一起而提高了可读性。
 *
 * 在初始化完成以后，由于给定的被执行函数的参数类型为void*，所以将参数传入为nullptr即可。
 * 整个过程中没有在实例化中给对象名字，所以要将两个括号连续使用，如果假设实例化以后对象的
 * 名字是tmp，那么该宏的工作过程实际上是如下两个语句：
 * Timer<void*, void*> tmp(void*, "tip");
 * tmp(nullptr);
 * 对比返回值不是void类型的函数，如果不用TIMING宏，那么计时功能的用法如下：
 * 假设需要被装饰的函数sum为：
 * int sum(int a, int b){
 *   return a + b;
 * }
 * 那么Timer类的用法如下：
 * Timer<int, int, int> timer(sum, "tip");
 * timer(1,2);
 * 可以看到，直接使用Timer类时被装饰的函数参数与调用反而被分开了，并且由于被装饰以后函数
 * 名的变化，远不如使用TIMING宏直观。
 */
#define TIMING(content, text)    \
    Timer<void *, void *>(       \
        [&](void *_) -> void * { \
            do {                 \
                content          \
            } while (false);     \
            return nullptr;      \
        },                       \
        text)(nullptr)

/**
 * @brief 统计函数运行时间的类
 *
 * @tparam RET 需要统计运行时间的函数返回值
 * @tparam ARGS 需要统计运行时间的函数参数
 */
template <typename RET, typename... ARGS>
class Timer : public Decorator<RET, ARGS...> {
   private:
    // 针对计时器的说明字符串（即在运行时间前针对该时间到底是哪个功能消耗时间的说明）
    std::string tip;

   public:
    explicit Timer(const std::function<RET(ARGS...)> &p_func,
                   std::string tip = "")
        : Decorator<RET, ARGS...>(p_func), tip(std::move(tip)) {}

    explicit Timer(RET (*p_func)(ARGS...), std::string tip = "")
        : Decorator<RET, ARGS...>(p_func), tip(std::move(tip)) {}

    RET operator()(ARGS... args) override {
        // 记录函数运行前的时间
        auto start = std::chrono::steady_clock::now();

        // 调用被装饰的函数
        RET result = this->decorated_func(args...);

        // 记录函数运行后的时间
        auto end = std::chrono::steady_clock::now();
        // 计算时间差并输出
        auto duration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start);

        std::cout << tip
                  << static_cast<double>(duration.count()) *
                         std::chrono::microseconds::period::num /
                         std::chrono::microseconds::period::den
                  << " s" << std::endl;

        return result;
    }
};
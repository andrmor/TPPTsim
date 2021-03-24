#ifndef out_h
#define out_h

#include <iostream>

template <typename Arg>
void outHelper(Arg&& arg)
{
    std::cout << arg << " ";
}

void out()
{
    std::cout << std::endl;
}

template <typename Arg, typename... Args>
void out(Arg&& arg, Args&&... rest)
{
    outHelper(arg);
    out(rest...);
}

#endif // out_h

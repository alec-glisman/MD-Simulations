//
// Created by Prakhar Srivastav (https://github.com/prakhar1989/progress-cpp)
// Output slightly modified to include more brackets and iterations / second by Alec Glisman on 1-5-12

#ifndef PROGRESSBAR_PROGRESSBAR_HPP
#define PROGRESSBAR_PROGRESSBAR_HPP

#include <chrono>
#include <iostream>

class ProgressBar {
private:
    unsigned int ticks = 0;

    const unsigned int total_ticks;
    const unsigned int bar_width;
    const char complete_char = '=';
    const char incomplete_char = ' ';
    const std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

public:
    ProgressBar(unsigned int total, unsigned int width, char complete, char incomplete) :
            total_ticks {total}, bar_width {width}, complete_char {complete}, incomplete_char {incomplete} {}

    ProgressBar(unsigned int total, unsigned int width) : total_ticks {total}, bar_width {width} {}

    unsigned int operator++() { return ++ticks; }

    void display() const
    {
        float progress = (float) ticks / total_ticks;
        auto pos = (int) (bar_width * progress);

        std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
        auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now-start_time).count();

        std::cout << "[";

        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) std::cout << complete_char;
            else if (i == pos) std::cout << ">";
            else std::cout << incomplete_char;
        }
        std::cout << "] " << int(progress * 100.0) << "% ["
                  << std::setprecision(2) << std::fixed << float(time_elapsed) / 1000.0 << "s, "
                  << std::setprecision(0) << std::fixed << float(ticks) / (float(time_elapsed) / 1000.0)
                  << "it/s"
                  << "]\r";
        std::cout.flush();
    }

    void done() const
    {
        display();
        std::cout << std::endl;
    }
};

#endif //PROGRESSBAR_PROGRESSBAR_HPP

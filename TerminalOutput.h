/*
 * The functions getTerminalWidth and printProgress
 * were taken from examples/md-flexible/src/Simulation.cpp
 */
#pragma once

#include <iostream>
#include <string>
#include <sys/ioctl.h>
#include <unistd.h>

class TerminalOutput {
private:
  size_t terminalWidth;
  bool showProgressBar;

public:
  TerminalOutput() {
    terminalWidth = getTerminalWidth();
    showProgressBar = true;
  }

  void setShowProgressBarValue(bool value) { showProgressBar = value; }

  size_t getTerminalWidth() {
    size_t terminalWidth = 0;
    // test all std pipes to get the current terminal width
    for (auto fd : {STDOUT_FILENO, STDIN_FILENO, STDERR_FILENO}) {
      if (isatty(fd)) {
        struct winsize w {};
        ioctl(fd, TIOCGWINSZ, &w);
        terminalWidth = w.ws_col;
        break;
      }
    }

    // if width is still zero try the environment variable COLUMNS
    if (terminalWidth == 0) {
      if (auto *terminalWidthCharArr = std::getenv("COLUMNS")) {
        // this pointer could be used to detect parsing errors via terminalWidthCharArr == end
        // but since we have a fallback further down we are ok if this fails silently.
        char *end{};
        terminalWidth = std::strtol(terminalWidthCharArr, &end, 10);
      }
    }

    // if all of the above fail fall back to a fixed width
    if (terminalWidth == 0) {
      // this seems to be the default width in most terminal windows
      terminalWidth = 80;
    }

    return terminalWidth;
  }
  void printProgress(size_t iterationProgress, size_t maxIterations) {
    if (not showProgressBar) { return; }

    // percentage of iterations complete
    const double fractionDone = static_cast<double>(iterationProgress) / static_cast<double>(maxIterations);

    // length of the number of maxIterations
    const int numCharsOfMaxIterations = static_cast<int>(std::to_string(maxIterations).size());

    // trailing information string
    std::stringstream info;
    info << std::setw(3) << std::round(fractionDone * 100) << "% " << std::setw(numCharsOfMaxIterations)
         << iterationProgress << "/";
    info << maxIterations;

    // actual progress bar
    std::stringstream progressbar;
    progressbar << "[";

    // the bar should fill the terminal window so subtract everything else (-2 for "] ")
    const int maxBarWidth = static_cast<int>(terminalWidth - info.str().size() - progressbar.str().size() - 2ul);
    // sanity check for underflow
    if (maxBarWidth > terminalWidth) {
      std::cerr << "Warning! Terminal width appears to be too small or could not be read. Disabling progress bar."
                << std::endl;
      setShowProgressBarValue(false);
      return;
    }
    const auto barWidth =
        std::max(std::min(static_cast<decltype(maxBarWidth)>(maxBarWidth * (fractionDone)), maxBarWidth), 1);
    // don't print arrow tip if >= 100%
    if (iterationProgress >= maxIterations) {
      progressbar << std::string(barWidth, '=');
    } else {
      progressbar << std::string(barWidth - 1, '=') << '>' << std::string(maxBarWidth - barWidth, ' ');
    }
    progressbar << "] ";
    // clear current line (=delete previous progress bar)
    std::cout << std::string(terminalWidth, '\r');
    // print everything
    std::cout << progressbar.str() << info.str() << std::flush;
  }
};

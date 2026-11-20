/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#include <cstdio>
#include <string>
#include <vector>
#include "Configuration.h"
#include "Testing.h"

Configuration parse(std::vector<std::string> args) {
    std::vector<char *> argv;
    for (std::string &arg : args) {
        argv.push_back(arg.data());
    }
    Configuration configuration;
    int error_code;
    ASSERT(configuration.parse(static_cast<int>(argv.size()), argv.data(), error_code));
    return configuration;
}

int main() {
    ASSERT(parse({"/opt/bin/empi", "--gabor", "in.bin", "out.db"}).variant == "SMP");
    ASSERT(parse({"empi", "--gabor", "--mmp3", "in.bin", "out.db"}).variant == "MMP3");
    ASSERT(parse({"empi", "--gabor", "in.bin", "out.db"}).signal_info.empty());
    ASSERT(parse({"empi", "--gabor", "--signal-info", "montage=CAR", "in.bin", "out.db"}).signal_info == "montage=CAR");

    ASSERT(parse({"empi", "--gabor", "in.bin", "out.db"}).decomp_argv.empty());
    ASSERT(parse({"/opt/bin/empi", "-i", "100", "--gabor", "--decomp-argv", "in.bin", "out.db"}).decomp_argv
           == "empi -i 100 --gabor --decomp-argv in.bin out.db");
    ASSERT(parse({"empi", "--gabor", "--decomp-argv", "--residual-log-dir=/home/user/logs/", "/home/user/data/my signal.bin", "C:\\books\\out.db"}).decomp_argv
           == "empi --gabor --decomp-argv --residual-log-dir=(...)/logs (...)/'my signal.bin' (...)/out.db");
    ASSERT(parse({"empi", "--gabor", "--decomp-argv", "--signal-info", "unit=fT, it's MEG", "in.bin", "out.db"}).decomp_argv
           == "empi --gabor --decomp-argv --signal-info 'unit=fT, it'\\''s MEG' in.bin out.db");
    ASSERT(parse({"empi", "--gabor", "--decomp-argv", "--residual-log-dir", "logs/", "/", "out.db"}).decomp_argv
           == "empi --gabor --decomp-argv --residual-log-dir logs/ / out.db");
    puts("OK");
}

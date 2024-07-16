//
// Created by Frank Yang on 7/16/24.
//

#ifndef CBU_CONSOLE_H
#define CBU_CONSOLE_H

#include <iostream>
#include "CBU_Balancer.h"

class CBU_Console {
public:
    CBU_Console() { std::cout << CBU_Console::version() << CBU_Balancer::version() << std::endl << CBU_Console::guide_message() << std::endl; }
    static std::string version();
    static std::string guide_message();
    void boot();
};

inline std::string CBU_Console::version() {
    std::string welcome_message = "";
    welcome_message += "Welcome to ChemicalBalancingUtilityNew Console v1.0.0!\n";
    welcome_message += "Programmed by Frank Yang in Jul, 2024 \n";
    return welcome_message;
}

inline std::string CBU_Console::guide_message() {
    std::string guide_message = "";
    guide_message += "Use quit() to quit the console.\n";
    guide_message += "Use multiple_results(off)` to disable multiple results, use multiple_results(on) to allow it. Multiple results is enabled by default.\n" ; 
    guide_message += "Directly type your chemical equation to call the built-in ChemicalBalancingUtility to balance.\n";
    return guide_message;
}

void CBU_Console::boot() {
    std::string command;
    CBU_Balancer balancer = CBU_Balancer();
    balancer.set_log_status(false);
    while (std::getline(std::cin, command)) {
        if (command == "quit()") {
            break;
        }
        else if (command == "multiple_results(on)") {
            balancer.set_multiple_results(true);
            std::cout << "Multiple results is allowed." << std::endl;
        }
        else if (command == "multiple_results(off)") {
            balancer.set_multiple_results(false);
            std::cout << "Multiple results is disabled. The balancing result my be not correct!" << std::endl;
        }
        else {
            balancer.balance(command);
            std::cout << balancer.get_result() << std::endl;
            balancer.clear_data();
        }
    }
}

#endif // CBU_CONSOLE_H

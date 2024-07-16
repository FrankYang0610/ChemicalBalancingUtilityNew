//
// Created and edited by Frank Yang on 7/14/24.
//

// Copyright (c) 2024 Frank Yang
// All rights reserved.
//
// This code is proprietary and confidential.
// Unauthorized copying of this file, via any medium is strictly prohibited.


// In this program (i.e., its context),
// A compound is a chemical compound, e.g., CaSO4, etc.
// A compound_str is a compound in std::string form, e.g., "CaSO4", etc.
// An entity is a chemical entity, including atom, ion and others, e.g., Ca, H, SO4, etc.
// An entity_str is an chemical entity (**with its coefficient**) in std::string form, e.g., "Ca", "H2", "SO4", etc.
// A compound composition is a std::map that maps the elements in the compound with their coefficients, e.g., "Ca"->1, "H"->2, "SO4"->1, etc.
// Use . to substitute · (in complex compounds e.g., CuSO4·5H2O)


#ifndef CBU_BALANCER_H
#define CBU_BALANCER_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <deque>
#include <set>
#define HAS_EXCEPTIONS
#define DEFAULT_MAX_COEF 20

class CBU_Balancer {
private:
    // Class settings
    bool _multiple_results; // will multiple results be allowed?
    unsigned _max_coef; // the maximum coefficient number allowed in the chemical equation
    bool _log_status; // Will there be loggs (only available when multiple_results is on)

    // Useful members
    std::pair<std::vector<std::string>, std::vector<std::string>> _reactants_and_products;
    ////
    //// These two following vectors should be used together
    std::vector<std::string> _elements; // Elements (in the recent equation)
    std::vector<std::vector<int>> _main_matrix; // Main matrix (reactants and products matrix, row = each element, column = each compound)
    ////
    std::vector<std::vector<unsigned>> _results_coefs;

    // Private methods
    static bool _is_valid_char(const char& c);
    static std::deque<std::string> _compound_str_separator(const std::string& compound_str);
    static std::pair<std::string, unsigned> _entity_str_to_entity_and_coef(const std::string& entity_str);
    static std::map<std::string, unsigned> _get_compound_composition(const std::string& compound_str);
    static std::vector<std::string> get_elements_from_compounds_composition(const std::vector<std::map<std::string, unsigned>>& compounds_composition);
    static std::vector<std::vector<unsigned>> _build_matrix(const std::vector<std::map<std::string, unsigned>>& compounds_composition, const std::vector<std::string>& elements);
    bool _solving_matrix_using_recursion(std::vector<unsigned>& coefficients_temporary, size_t floor);

    //// Math tools
    static bool _are_linear_dependent(const std::vector<unsigned>& a, const std::vector<unsigned>& b);
    void _filter_linear_independent_results();

    //// String tools
    static std::string _remove_spaces(const std::string& str);
    static std::pair<std::string, std::string> _get_reactants_and_products_str(const std::string& equation);
    static std::vector<std::string> _separate_half_equation(const std::string& half_equation);
    static std::pair<std::vector<std::string>,std::vector<std::string>> _get_compounds_str(const std::string& equation);
    void _balance_with_given_compounds_str(const std::vector<std::string>& reactants, const std::vector<std::string>& products);
public:
    // Constructors, getters and setters
    CBU_Balancer() : _multiple_results(true), _max_coef(DEFAULT_MAX_COEF), _log_status(true), _reactants_and_products(),_elements(),_main_matrix(), _results_coefs() { };
    void set_multiple_results(bool option) { this->_multiple_results = option; }
    void set_max_coef(unsigned max_coef) { this->_max_coef = max_coef;  }
    void set_log_status(bool option) { this->_log_status = option; };
    // Public interfaces
    void balance_with_given_compounds(const std::vector<std::string>& reactants, const std::vector<std::string>& products) { // alias
        _balance_with_given_compounds_str(reactants, products);
    }
    void balance(const std::string& equation);
    // Common getters
    std::pair<std::vector<std::string>, std::vector<std::string>> get_reactants_and_products() { return this->_reactants_and_products; }
    std::pair<std::vector<std::string>, std::vector<std::vector<int>>> get_main_matrix() { return {this->_elements, this->_main_matrix}; };
    [[nodiscard]] std::string get_result() const;
    void clear_data();
    static std::string version();
};

// This method checks if a character (const char &c) is valid in a chemical equation
inline bool CBU_Balancer::_is_valid_char(const char& c) {
    return std::isupper(c) || std::islower(c) || std::isdigit(c) || (c == '(') || (c == ')');
}

// This method takes a compound string (const std::string &compound_str) as input and returns the chemical entity strings (i.e., chemical entity + coefficients strings) in the compound.
HAS_EXCEPTIONS
std::deque<std::string> CBU_Balancer::_compound_str_separator(const std::string& compound_str) {
    try {
        std::deque<std::string> compound_separated; // Note: CaSO4 will be separated into ["Ca", "SO4"] and stored here
        int in_parentheses = 0;

        for (const auto& c : compound_str) {
            if (!_is_valid_char(c)) { throw std::runtime_error("INVALID_CHAR"); }

            if (c == '(') {
                in_parentheses++;
                if (in_parentheses == 1) { compound_separated.emplace_back(1, c); continue; } // 'actually' not in parentheses
            }

            if (in_parentheses == 0) {
                if (std::isupper(c)) { compound_separated.emplace_back(1, c); continue; } // compound_separated.emplace_back(std::string(1, c));

                // c is not a capital letter.
                if (!compound_separated.empty()) {
                    if (std::islower(c) || std::isdigit(c)) { compound_separated.back() += c; }
                    else { throw std::runtime_error(compound_str + ": " + "UNMATCHED_PARENTHESES"); }
                } else { throw std::runtime_error("INVALID_CHAR_OR_CHAR_POSITION");  } // for example, compound_str == "c26" or ")CO2"
            }

            else if (in_parentheses > 0) { // in_parentheses
                if (c == ')') { in_parentheses--; }

                if (!compound_separated.empty()) { compound_separated.back() += std::string(1, c); }
                else { throw std::runtime_error("EMPTY_CONTAINER"); }
            }

            else { throw std::runtime_error("UNMATCHED_PARENTHESES"); } // in_parentheses < 0
        }

        if (in_parentheses) { throw std::runtime_error("UNMATCHED_PARENTHESES"); }

        return compound_separated;
    } catch (std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
}

// This method converts a chemical entity string into an std::pair of <chemical entity, coefficient>.
HAS_EXCEPTIONS
inline std::pair<std::string, unsigned> CBU_Balancer::_entity_str_to_entity_and_coef(const std::string& entity_str) {
    try {
        if (entity_str.empty()) {
            throw std::runtime_error("EMPTY_ENTITY_STR");
        }

        bool is_subcompound = (entity_str[0] == '(');
        int separating_index = -1; // the beginning position index of the coefficient

        if (std::isdigit(entity_str[entity_str.size() - 1])) {
            for (size_t index = entity_str.size() - 1; index >= 0; index--) {
                if (!isdigit(entity_str[index])) {
                    separating_index = static_cast<int>(index) + 1;
                    break;
                }
            }
        }

        std::pair <std::string, unsigned> entity_and_coefficient;

        // To enhance readability, these conditional statements will not be shortened by the ternary operator (?:).
        if (!is_subcompound) {
            if (separating_index == -1) { // coef = 1
                entity_and_coefficient.first = entity_str;
                entity_and_coefficient.second = 1;
            } else {
                entity_and_coefficient.first = entity_str.substr(0, separating_index);
                entity_and_coefficient.second = std::stoi(entity_str.substr(separating_index, entity_str.size() - separating_index));
            }
        } else { // is_subcompound (within parentheses)
            if (separating_index == -1) { // coef = 1
                entity_and_coefficient.first = entity_str.substr(1, entity_str.size() - 2);
                entity_and_coefficient.second = 1;
            } else {
                entity_and_coefficient.first = entity_str.substr(1, separating_index - 2);
                entity_and_coefficient.second = std::stoi(entity_str.substr(separating_index, entity_str.size() - separating_index));
            }
        }

        if (!is_subcompound && entity_and_coefficient.first.size() >= 3) {
            std::cout << ("WARNING: special element \"" + entity_and_coefficient.first + "\", the equation may not exist.") << std::endl;
        }

        return entity_and_coefficient;
    } catch (const std::invalid_argument &ia) {
        std::cout << ia.what() << std::endl;
        return {};
    } catch (const std::out_of_range &oor) {
        std::cout << oor.what() << std::endl;
        return {};
    } catch (const std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
} // .first is the entity, .second is its coefficient

// This method takes a compound string as input and outputs a std::map that records the coefficient of each element in the compound.
// The compound string must be valid, or there will be exceptions.
HAS_EXCEPTIONS
std::map<std::string, unsigned> CBU_Balancer::_get_compound_composition(const std::string& compound_str) {
    try {
        std::deque<std::string> entities_str = CBU_Balancer::_compound_str_separator(compound_str); // Note: CaSO4 will be separated into ["Ca", "SO4"] and stored here
        if (entities_str.empty()) { throw std::runtime_error("INVALID_COMPOUND"); }

        std::map<std::string, unsigned> composition; // elements and their counts in [current] compound

        for (const auto& entity_str : entities_str) {
            std::pair<std::string, int> entity_and_coefficient = CBU_Balancer::_entity_str_to_entity_and_coef(entity_str);

            if (entity_str[0] == '(') {
                std::map <std::string, unsigned> composition_subcompound = CBU_Balancer::_get_compound_composition(entity_and_coefficient.first);
                for (const auto &cs : composition_subcompound) {
                    composition[cs.first] += cs.second * entity_and_coefficient.second;
                }
            } else {
                composition[entity_and_coefficient.first] += entity_and_coefficient.second;
            }
        }
        return composition;
    } catch (const std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
}

// This method takes a compound composition (std::map) vector (i.e., std::vector<std::map>) as input and outputs the elements list (sorted) in these compounds
HAS_EXCEPTIONS
std::vector<std::string>
CBU_Balancer::get_elements_from_compounds_composition(const std::vector<std::map<std::string, unsigned>>& compounds_composition) {
    try {
        if (compounds_composition.empty()) {
            throw std::runtime_error("EMPTY_CONTAINER");
        }

        std::set<std::string> elements_set;
        for (const auto& compound_composition : compounds_composition) {
            for (const auto& element : compound_composition) {
                elements_set.insert(element.first);
            }
        }

        std::vector<std::string> elements = {elements_set.begin(), elements_set.end()};
        std::sort(elements.begin(), elements.end());
        return elements;
    } catch (std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
}

// This method takes compounds composition (std::vector<std::map>) and elements list (std::vector) as input and outputs a matrix, where [the row is each element] and [the column is each compound].
HAS_EXCEPTIONS
std::vector<std::vector<unsigned>>
CBU_Balancer::_build_matrix(const std::vector<std::map<std::string, unsigned>> &compounds_composition, const std::vector<std::string>& elements) {
    try {
        if (compounds_composition.empty() || elements.empty()) {
            throw std::runtime_error("EMPTY_CONTAINER");
        }

        // 2D vector
        std::vector<std::vector<unsigned>> matrix(elements.size(), std::vector<unsigned>(compounds_composition.size(), 0));

        for (size_t i = 0; i < elements.size(); i++) { // row: element
            for (size_t j = 0; j < compounds_composition.size(); j++) { // column: compounds
                matrix[i][j] = static_cast<int>(compounds_composition[j].count(elements[i]) ? compounds_composition[j].at(elements[i]) : 0);
            }
        }

        return matrix;
    } catch (std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
}

// This method is used to solve the matrix (where, the result is a std::vector<unsigned>) using recursion strategy. To the inputs, coefficients is a temporary vector and floor is the recursion floor.
bool CBU_Balancer::_solving_matrix_using_recursion(std::vector<unsigned>& coefficients_temporary, size_t floor) {
    if (floor == this->_main_matrix[0].size()) { // recursion end
        bool balanced = true;
        for (const auto& element :this->_main_matrix) {
            int sum = 0;
            for (size_t i = 0; i < coefficients_temporary.size(); i++) {
                sum += static_cast<int>(coefficients_temporary[i]) * element[i];
            }
            if (sum != 0) { balanced = false; break; }
        }

        if (balanced) {
            if (this->_multiple_results && this->_log_status) {
                std::clog << "A possible result found. " << std::endl;
            }
            if (this->_multiple_results || this->_results_coefs.empty()) {
                this->_results_coefs.push_back(coefficients_temporary);
            }
            return true;
        }
        return false;
    }

    for (size_t i = 1; i <= this->_max_coef; i++) {
        coefficients_temporary[floor] = i;
        bool balanced = CBU_Balancer::_solving_matrix_using_recursion(coefficients_temporary, floor + 1);
        if (balanced && !this->_multiple_results) { return true; }
    }

    return false;
}

// This method tests if two vectors are linear [dependent].
inline bool CBU_Balancer::_are_linear_dependent(const std::vector<unsigned> &a, const std::vector<unsigned> &b) {
    if (a.empty() || b.empty() || (a.size() != b.size())) { return false; }
    double ratio = 0;
    bool ratio_set = false;

    for (size_t i = 0; i < a.size(); i++) {
        if (a[i] != 0) {
            double current_ratio = static_cast<double>(b[i]) / static_cast<double>(a[i]);
            if (!ratio_set) {
                ratio = current_ratio;
                ratio_set = true;
            } else if (std::abs(ratio - current_ratio) > 1e-9) {
                return false;
            }
        } else if (b[i] != 0) { // a[i] == 0 && b[i] != 0
            return false;
        }
    }
    return true;
}

// This method removes linear dependent items.
// Modifying private members
void CBU_Balancer::_filter_linear_independent_results() {
    if (this->_multiple_results && this->_log_status) { std::clog << "Filtering linear independent items" << std::endl; }

    std::vector<std::vector<unsigned>> independent_results_coefs;
    for (const auto& result : this->_results_coefs) {
        bool independent = true;
        if (!independent_results_coefs.empty()) {
            for (const auto& confirmed_result : independent_results_coefs) {
                if (CBU_Balancer::_are_linear_dependent(confirmed_result, result)) {
                    independent = false;
                }
            }
        }
        if (independent) {
            independent_results_coefs.push_back(result);
        }
    }
    this->_results_coefs = std::move(independent_results_coefs);
}

// This method removes all spaces from a string and returns the new string.
inline std::string CBU_Balancer::_remove_spaces(const std::string &str) {
    std::string result = str;
    result.erase(std::remove_if(result.begin(), result.end(), [](unsigned char c) { return std::isspace(c); }), result.end());
    return result;
}

// This method generates reactants and products string (std::pair<std::string, std::string> from an equation string)
HAS_EXCEPTIONS
std::pair<std::string, std::string> CBU_Balancer::_get_reactants_and_products_str(const std::string &equation) {
    try {
        std::string delimiter = "->";
        std::vector<std::string> tokens;
        size_t start = 0;
        size_t end = equation.find(delimiter);

        while (end != std::string::npos) {
            tokens.push_back(equation.substr(start, end - start));
            start = end + delimiter.length();
            end = equation.find(delimiter, start);
        }
        tokens.push_back(equation.substr(start));

        if (tokens.size() == 2) {
            return {tokens[0], tokens[1]};
        } else {
            throw std::runtime_error("INVALID_EQUATION");
        }
    } catch (const std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
}

// This method takes a half equation (for example, the half equations of "N + O2 -> NO2" is "N + O2" and "NO2") as input and outputs the compounds str of that half equation.
// Generally, 'half_equations' are from the previous method, _get_reactants_and_products_str.
// View _get_compounds_str
HAS_EXCEPTIONS
std::vector<std::string> CBU_Balancer::_separate_half_equation(const std::string &half_equation) {
    try {
        if (half_equation.empty()) {
            throw std::runtime_error("INCOMPLETE_EQUATION");
        }
        std::string delimiter = "+";
        std::vector<std::string> tokens;
        size_t start = 0;
        size_t end = half_equation.find(delimiter);

        while (end != std::string::npos) {
            tokens.push_back(half_equation.substr(start, end - start));
            start = end + delimiter.length();
            end = half_equation.find(delimiter, start);
        }
        tokens.push_back(half_equation.substr(start));
        return tokens;
    } catch (const std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
}

// This method is used to generate the compounds str (std::pair for reactant and products; std::vector<std::string> for compounds str) from an equation
HAS_EXCEPTIONS
std::pair<std::vector<std::string>, std::vector<std::string>>
CBU_Balancer::_get_compounds_str(const std::string &equation) {
    try {
        std::string equation_without_space = CBU_Balancer::_remove_spaces(equation);
        std::pair<std::string, std::string> reactants_and_products_str = CBU_Balancer::_get_reactants_and_products_str(equation_without_space);

        std::pair<std::vector<std::string>, std::vector<std::string>> compounds_str = {
            CBU_Balancer::_separate_half_equation(reactants_and_products_str.first),
            CBU_Balancer::_separate_half_equation(reactants_and_products_str.second)
        };
        return compounds_str;
    } catch (const std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
        return {};
    }
}

// This method is used to generate balancing result from given compounds
// Processing of complex compounds is included
HAS_EXCEPTIONS
void CBU_Balancer::_balance_with_given_compounds_str(const std::vector<std::string>& reactants, const std::vector<std::string>& products) {
    try {
        if (reactants.empty() || products.empty()) { throw std::runtime_error("EMPTY_COMPOUND_LIST"); }

        this->_reactants_and_products = {reactants, products};

        // Initialization && parse
        std::vector<std::map<std::string, unsigned>> reactants_composition, products_composition;
        for (const auto& reactant : reactants) {
            if (reactant.empty()) { throw std::runtime_error("INCOMPLETE_REACTANT"); }
            reactants_composition.push_back(CBU_Balancer::_get_compound_composition(reactant));
        } // to all reactants
        for (const auto& product : products) {
            if (product.empty()) { throw std::runtime_error("INCOMPLETE_REACTANT"); }
            products_composition.push_back(CBU_Balancer::_get_compound_composition(product));
        } // to all products

        std::vector<std::string> reactant_elements = CBU_Balancer::get_elements_from_compounds_composition(reactants_composition);
        std::vector<std::string> product_elements = CBU_Balancer::get_elements_from_compounds_composition(products_composition);

        if (reactant_elements.size() == product_elements.size()) {
            auto mismatch_result = std::mismatch(reactant_elements.begin(), reactant_elements.end(),
                                                 product_elements.begin(), product_elements.end());
            if (mismatch_result.first != reactant_elements.end() || mismatch_result.second != product_elements.end()) {
                throw std::runtime_error("ELEMENTS_MISMATCH_BETWEEN_REACTANTS_AND_PRODUCTS");
            }
        } else {
            throw std::runtime_error("ELEMENTS_MISMATCH_BETWEEN_REACTANTS_AND_PRODUCTS");
        }

        this->_elements = std::move(reactant_elements); // save to private member

        // Build matrices
        std::vector<std::vector<unsigned>> reactant_matrix = CBU_Balancer::_build_matrix(reactants_composition, _elements);
        std::vector<std::vector<unsigned>> product_matrix = CBU_Balancer::_build_matrix(products_composition, _elements);

        std::vector<std::vector<int>> main_matrix(
            this->_elements.size(), // row is the element
            std::vector<int>(reactants.size() + products.size(), 0) // column is the reactants and the products
            // all filling 0
        );

        for (size_t i = 0; i < this->_elements.size(); i++) {
            for (size_t j = 0; j < reactants.size(); ++j) {
                main_matrix[i][j] = static_cast<int>(reactant_matrix[i][j]);
            }
            for (size_t j = 0; j < products.size(); j++) {
                main_matrix[i][reactants.size() + j] = - static_cast<int>(product_matrix[i][j]);
            }
        }

        this->_main_matrix = std::move(main_matrix); // save to private member

        // Brute force balancing
        std::vector<unsigned> results_coefficients(reactants.size() + products.size(), 0);
        CBU_Balancer::_solving_matrix_using_recursion(results_coefficients, 0); // !
        CBU_Balancer::_filter_linear_independent_results(); // here will be faster

        if (this->_results_coefs.empty()) {
            throw std::runtime_error("FAILED_TO_BALANCE");
        }
    } catch (std::runtime_error &re) {
        std::cerr << re.what() << std::endl;
    }
}

// Main balance method
void CBU_Balancer::balance(const std::string &equation) {
    std::pair<std::vector<std::string>, std::vector<std::string>> reactants_and_products = _get_compounds_str(equation);
    this->_reactants_and_products = reactants_and_products;
    _balance_with_given_compounds_str(reactants_and_products.first, reactants_and_products.second);
}

// Get results
HAS_EXCEPTIONS
inline std::string CBU_Balancer::get_result() const {
    try {
        std::string result_str = "";
        if (this->_reactants_and_products.first.empty() || this->_reactants_and_products.second.empty()) {
            throw std::runtime_error("NO_RESULT");
        }

        unsigned reactants_count = _reactants_and_products.first.size();
        unsigned products_count = _reactants_and_products.second.size();

        if (this->_results_coefs.size() >= 2) { std::cout << "There are multiple possible solutions. Please select the most correct one." << std::endl; }

        for (const auto& solution : this->_results_coefs) { // solution: const std::vector<unsigned>
            for (size_t i = 0; i < reactants_count; i++) {
                auto coefficient = solution[i];
                result_str += ((coefficient != 1) ? std::to_string(coefficient) : "");
                result_str += _reactants_and_products.first[i];
                if (i < reactants_count - 1) {
                    result_str += " + ";
                }
            }
            result_str += " == ";
            for (size_t i = 0; i < products_count; i++) {
                auto coefficient = static_cast<unsigned>(solution[reactants_count + i]);
                result_str += ((coefficient != 1) ? std::to_string(coefficient) : "");
                result_str += _reactants_and_products.second[i];
                if (i < products_count - 1) {
                    result_str += " + ";
                }
            }
            result_str += '\n';
        }
        result_str = result_str.substr(0, result_str.size() - 1);
        return result_str;
    } catch (const std::runtime_error& re) {
        std::cerr << re.what() << std::endl;
        return "";
    }
}

inline void CBU_Balancer::clear_data() {
    this->_reactants_and_products.first.clear();
    this->_reactants_and_products.second.clear();
    this->_elements.clear();
    this->_main_matrix.clear();
    this->_results_coefs.clear();
}

// Get program version
inline std::string CBU_Balancer::version() {
    std::string version_str = "";
    version_str += "ChemicalBalancingUtility, Version New.1.0.0\n";
    version_str += "Programmed by Frank Yang in Jul 2024\n";
    return version_str;
}

#endif // CBU_BALANCER_H

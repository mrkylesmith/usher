/**
 * Copyright (c) 2006-09 Maciej F. Boni.  All Rights Reserved.
 * 
 * FILE NAME:   String.cpp
 * CREATED ON:  16 June 2011, 11:33
 * AUTHOR:      Ha Minh Lam
 * 
 * DESCRIPTION: The implementation of all the utilities for string processing.
 * 
 * HISTORY:     Version    Date         Description
 *              1.0        2011-06-16   created
 * 
 * VERSION:     1.0
 * LAST EDIT:   16 June 2011
 */

#include "String.h"


const bool String::isInteger(std::string str) {
    for (unsigned long i = 0; i < str.length(); i++) {
        if (!isdigit(str[i])) {
            return false;
        }
    }
    return true;
}

std::string String::formatInt(const long& integer, int minDigitNum) {
    std::stringstream stream;
    std::string result;

    stream << integer;
    result = stream.str();
    minDigitNum = minDigitNum - static_cast<int>(result.length());

    while (minDigitNum > 0) {
        result = "0" + result;
        minDigitNum--;
    }

    return result;
}

std::string String::trim(std::string inStr) {
    std::string::iterator it = inStr.begin();
    while (it != inStr.end()) {
        if (isspace(*it)) {
            inStr.erase(it);
        } else {
            break;
        }
    }

    it = inStr.end();
    it--;
    while (it != inStr.begin()) {
        if (isspace(*it)) {
            inStr.erase(it);
            it--;
        } else {
            break;
        }
    }

    return inStr;
}

std::string String::deleteAllSpace(std::string inStr) {
    std::string::iterator it = inStr.begin();
    while (it != inStr.end()) {
        if (isspace(*it)) {
            inStr.erase(it);
        } else {
            it++;
        }
    }

    return inStr;
}

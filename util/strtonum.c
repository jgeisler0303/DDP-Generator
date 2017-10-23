#include <cmath>
#include "strtonum.h"

int strtoint(const char* str, long *ret, const char **ret_ep) {
    if(ret==NULL || str==NULL) return 1;

    if(ret_ep) *ret_ep= str;
    if(str[0]=='\0') return 2;
    
    errno = 0;
    char *endptr;
    long val = strtol(str, &endptr, 10);
    *ret= val;
    if(ret_ep) *ret_ep= endptr;

    if (errno == ERANGE) {
        switch(val) {
        case LONG_MIN:
            return 3;
        case LONG_MAX:
            return 4;
        default:
            return 5;
        }
    } else if (errno != 0) {
        return 5;
//     } else if (*endptr != '\0') {
//         return 6;
    } else if (endptr == str) {
        for(; *str != '\0'; str++)
            if(!isspace(*str))
                return 7;
        return 8;
    }

    if(ret_ep) {
        while(isspace(*endptr) && *endptr!='\0') endptr++;
        *ret_ep= endptr;
    }
    
    return 0;
}

int strtodouble(const char* str, double *ret, const char **ret_ep) {
    if(ret==NULL || str==NULL) return 1;

    if(ret_ep) *ret_ep= str;
    if(str[0]=='\0') return 2;
    
    errno = 0;
    char *endptr;
    double val = strtod(str, &endptr);
    *ret= val;
    if(ret_ep) *ret_ep= endptr;

    if (errno == ERANGE) {
        if(val==HUGE_VAL)
            return 3;
        else
            return 5;
    } else if (errno != 0) {
        return 5;
//     } else if (*endptr != '\0') {
//         return 6;
    } else if (endptr == str) {
        for(; *str != '\0'; str++)
            if(!isspace(*str))
                return 7;
        return 8;
    }

    if(ret_ep) {
        while(isspace(*endptr) && *endptr!='\0') endptr++;
        *ret_ep= endptr;
    }
    
    return 0;
}

const char *strtonum_error(int i) {
    static const char strtoint_error1[]= "str or ret == NULL";
    static const char strtoint_error2[]= "empty string";
    static const char strtoint_error3[]= "underflow";
    static const char strtoint_error4[]= "overflow";
    static const char strtoint_error5[]= "unknown error";
    static const char strtoint_error6[]= "not entire string";
    static const char strtoint_error7[]= "no number found";
    static const char strtoint_error8[]= "all white space";
    
    switch(i) {
        case 1: return strtoint_error1;
        case 2: return strtoint_error2;
        case 3: return strtoint_error3;
        case 4: return strtoint_error4;
        case 5: return strtoint_error5;
        case 6: return strtoint_error6;
        case 7: return strtoint_error7;
        case 8: return strtoint_error8;
        default: return strtoint_error5;
    }
}
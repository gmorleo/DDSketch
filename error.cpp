//
// Created by giuseppe on 25/11/19.
//

#include "error.h"

int printError(int error) {

    switch (error) {
        case -1:
            cout << "Generic error" << endl;
            break;
        case -2:
            cout << "Not enough memory" << endl;
            break;
        case -3:
            cout << "Error while opening file" << endl;
            break;
        case -4:
            cout << "The sketch data structure is not valid" << endl;
    }

    return 0;
}
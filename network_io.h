#ifndef NETWORK_IO_H
#define NETWORK_IO_H

#include <iostream>

#include "network.h"


void read_network(std::istream & in, AbstractNetwork & net);
void write_network(std::ostream & on, const AbstractNetwork & net);




#endif	// NETWORK_IO_H

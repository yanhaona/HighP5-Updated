/*
 * cpuCores.h
 *
 *  Created on: Aug, 2024
 *      Author: yan
 */

#ifndef CPUCORES_H_
#define CPUCORES_H_

int* getCoreNumberingArray(const char *machine, int threadCount);
int getCoreCount(const char *machine, int threadCount);

#endif /* CPUCORES_H_ */

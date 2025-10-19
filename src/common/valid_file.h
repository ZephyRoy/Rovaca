#ifndef FILE_SYSTEM_MACRO_H
#define FILE_SYSTEM_MACRO_H

#include <fcntl.h>
#include <unistd.h>

#include <iostream>

bool isFileExists(const std::string& path) { return access(path.c_str(), F_OK) == 0; }
bool isFileReadable(const std::string& path) { return access(path.c_str(), R_OK) == 0; }
bool isFileWritable(const std::string& path) { return access(path.c_str(), W_OK) == 0; }

#define CHECK_FILE_EXIST(path)                             \
    do {                                                   \
        if (!isFileExists(path)) {                         \
            RovacaLogger::error("file: {} not exist", path); \
            return false;                                  \
        }                                                  \
    } while (0)

#define CHECK_FILE_READALBE(path)                             \
    do {                                                      \
        if (!isFileReadable(path)) {                          \
            RovacaLogger::error("file: {} not readable", path); \
            return false;                                     \
        }                                                     \
    } while (0)

#define CHECK_FILE_WRITEABLE(path)                             \
    do {                                                       \
        int fd = open(path, O_WRONLY | O_CREAT, 0666);         \
        if (fd == -1) {                                        \
            RovacaLogger::error("file: {} not writeable", path); \
            return false;                                      \
        }                                                      \
        close(fd);                                             \
    } while (0)

#endif  // FILE_SYSTEM_MACRO_H
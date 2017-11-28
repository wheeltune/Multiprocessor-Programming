#ifndef LOG_H
#define LOG_H

static int ENABLE_LOGGING = 0;

void message(const char *text) {
    if (ENABLE_LOGGING == 1) {
        printf("%s\n", text);
        fflush(stdout);
    }
}

void messageD(const char *text, const char *detail) {
    if (ENABLE_LOGGING) {
        printf("%s: %s\n", text, detail);
        fflush(stdout);
    }
}

#endif // LOG_H

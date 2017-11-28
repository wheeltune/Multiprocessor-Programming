#ifndef LOG_H
#define LOG_H

static int ENABLE_LOGGING = 0;

void message(const char *text);
void messageD(const char *text, const char *detail);

#endif // LOG_H

CC ?= cc
CFLAGS ?= -Wall
LIBS = -lpng -lm
LDFLAGS ?= -s
RM = rm -f
TARGET = xpng

all: $(TARGET)

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LIBS) -o $@

$(TARGET).o: $(TARGET).c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(TARGET).o

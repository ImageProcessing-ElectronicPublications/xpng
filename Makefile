CC ?= cc
CFLAGS ?= -Wall
LIBS = -lpng -lz -lm
LDFLAGS ?= -s
RM = rm -f
TARGET = xpng

all: $(TARGET)

$(TARGET): src/$(TARGET).o
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(LIBS) -o $@

src/$(TARGET).o: src/$(TARGET).c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) src/$(TARGET).o

MAIN := main.pdf
CS := cheatsheet.pdf
CHAPS := $(wildcard chap*.typ)

all: $(MAIN) $(CS)

$(MAIN): main.typ $(CHAPS)
	typst compile main.typ

$(CS): cheatsheet.typ
	typst compile cheatsheet.typ

clean:
	rm -f $(MAIN) $(CS)

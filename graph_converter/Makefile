all: graph_formater

SRC = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SRC))


#OBJECTS = [fasta_util.o message_value.o \
	shingling_cluster_generator_optimized.o nmi_messages.o \
	label_as_value.o coupled_labels.o matched_cluster.o mr_util.o \
	/home/armen.abnousi/shingle/merge_evaluator/merge_evaluator.o \
	union_find.o mpi_util.o fvalue_evaluator.o graph_formater.o shingling_toolbox.o]

#%.o: %.cpp
#	echo "executing o"; \
#	echo $(shell pwd); \
#	$(CPP) -c $(CPPFLAGS)  $<

graph_formater: $(OBJECTS)
	$(LINK) $(LINKFLAGS) $(OBJECTS) ../src/graph_formater.o ../src/fasta_util.o ../src/mr_util.o $(USRLIB) $(SYSLIB) -o ../bin/graph_formatter

#%.o: %.cpp
#	$(CPP) -c -o $@ $< $(CPPFLAGS)

#%.o: %.cpp
#	echo "executing o"; \
#	echo $(pwd); \
#	$(CPP) $(CPPFLAGS) -c $< -I$(INCLUDES) -o $@


#%.o:
#	$(CC) $(CFLAGS) -c $< -I$(INCLUDES) -o $@

#%.o:%.cpp
#	cd src/; \
#	$(CPP) $(CPPFLAGS) -c $< -I$(INCLUDES) -o $@

%.o: %.cpp
	$(CPP) -c $(CPPFLAGS) $<

clean:
	rm -f $(TARGET).o $(OBJECTS)

wipe:
	rm -f $(TARGET).o $(OBJECTS) $(TARGET) *~ *.bak

message:
	echo "Executables: " $(TARGET) " have been created"

CXX=g++
CXXFLAGS=-O2 -g -march=native -mtune=native
LIBS=-lduneanaobj_StandardRecordProxy -lduneanaobj_StandardRecord -lduneanaobj_StandardRecord_dict
INC_DIR=-I${SRPROXY_INC} -I${DUNEANAOBJ_INC}
LIB_DIR=-L${DUNEANAOBJ_LIB}
ROOT_CFLAGS=`root-config --cflags`
ROOT_LIBS=`root-config --libs`

selection: purity_selection.cxx
	$(CXX) $(CXXFLAGS) $(INC_DIR) $(LIB_DIR) $(LIBS) $(ROOT_CFLAGS) $(ROOT_LIBS) $^ -o $@

clean:
	rm selection

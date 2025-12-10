// SWIG file CloudMesher.i

%{
#include "otmeshing/CloudMesher.hxx"
%}

%include CloudMesher_doc.i

%copyctor OTMESHING::CloudMesher;

%include otmeshing/CloudMesher.hxx

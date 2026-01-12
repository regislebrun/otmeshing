// SWIG file ConvexHullMesher.i

%{
#include "otmeshing/ConvexHullMesher.hxx"
%}

%include ConvexHullMesher_doc.i

%copyctor OTMESHING::ConvexHullMesher;

%include otmeshing/ConvexHullMesher.hxx

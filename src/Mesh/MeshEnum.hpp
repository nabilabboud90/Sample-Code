#ifndef MESHENUM_HPP_
#define MESHENUM_HPP_

namespace NSMesh {

enum Topology { tp_ELEMENT = 0, tp_FACE = 1, tp_EDGE = 2, tp_VERTEX = 3, tp_FIRST = tp_ELEMENT, tp_LAST = tp_VERTEX };
enum Ownership { o_OWNED = 0, o_SHARED = 1, o_UNIVERSAL = 2, o_FIRST = o_OWNED, o_LAST = o_UNIVERSAL };
enum class ElementTopology { et_UNDEFINED, et_LINE2, et_LINE3, et_TRI3, et_TRI6, et_TET4, et_TET10, et_QUAD4, et_QUAD8, et_HEX8, et_HEX20, et_FIRST = et_UNDEFINED, et_LAST = et_HEX20 };
enum class MeshQueryType { mqt_R0, mqt_R2, mqt_FIRST = mqt_R0, mqt_LAST = mqt_R2 };
  
} // namespace NSMesh

#endif

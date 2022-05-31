
#include "graph/abstract_edge.h"

namespace scnu
{
    /**
     * @details construct an edge with given source vertex id and destination vertex id
     * @param other_source_vertex_id
     * @param other_destination_vertex_id
     */
    abstract_edge::abstract_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id)
    :source_vertex_id(other_source_vertex_id),destination_vertex_id(other_destination_vertex_id){

    }

    /**
     * @details construct an edge with the given edge
     * @param other_edge
     */
    abstract_edge::abstract_edge(const std::shared_ptr<abstract_edge>& other_edge):
            abstract_edge(other_edge->get_source_vertex_id(),other_edge->get_destination_vertex_id()){
    }


    /**
     * @details get the source vertex id of this edge
     * @return
     */
    uint32_t abstract_edge::get_source_vertex_id() const {
        return this->source_vertex_id;
    }

    /**
     * @details get the destination vertex id of this vertex
     * @return
     */
    uint32_t abstract_edge::get_destination_vertex_id() const {
        return this->destination_vertex_id;
    }

    /**
     * @details swap the source vertex id and destination vertex id
     */
    void abstract_edge::swap(uint32_t source_id,uint32_t destination_id)
    {
        source_vertex_id = source_id;
        destination_vertex_id = destination_id;
    }

    /**
     * @details overload operator == to remove repeated edges
     * @param other_edge
     * @return
     */
    bool abstract_edge::operator==(const abstract_edge &other_edge) const{
        return source_vertex_id == other_edge.get_source_vertex_id()
               && destination_vertex_id == other_edge.get_destination_vertex_id();
    }
}




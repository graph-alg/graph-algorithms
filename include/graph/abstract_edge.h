/******************************************************************************************************************
* @copyright  : Boost Software License
* @file       : abstract_edge.h
* @brief      : An abstract edge structure
* @version    : 1.1
* @date       : 2020/10/15
******************************************************************************************************************/
#pragma once
#include "graph/graph_utility.h"

namespace scnu {
    /**
     * @details an abstract edge without direction and weight
     * @remarks source vertex and destination vertex have the sub type
     */
    class abstract_edge {
    public:
        explicit abstract_edge(uint32_t other_source_vertex_id, uint32_t other_destination_vertex_id);

        explicit abstract_edge(const std::shared_ptr<abstract_edge>& other_edge);

        virtual ~abstract_edge() = default;

        [[nodiscard]] uint32_t get_source_vertex_id() const;

        [[nodiscard]] uint32_t get_destination_vertex_id() const;

        void swap(uint32_t source_id,uint32_t destination_id);

        virtual bool operator == (const abstract_edge &other_edge) const;

    private:
        int source_vertex_id;
        int destination_vertex_id;
    };

    struct hash_abstract_edge {
        size_t operator()(const shared_ptr<abstract_edge> &e) const {
            stringstream input_stream;
            input_stream<<e->get_source_vertex_id()<< "," <<e->get_destination_vertex_id();
            return hash<std::string>()(input_stream.str());
        }
    };

    struct equal_abstract_edge {
        bool operator()(const shared_ptr<abstract_edge> &e1, const shared_ptr<abstract_edge> &e2) const {
            return e1->get_source_vertex_id() == e2->get_source_vertex_id()
                   && e1->get_destination_vertex_id() == e2->get_destination_vertex_id();
        }
    };
}



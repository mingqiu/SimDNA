//
// Created by Mingqiu Wang on 1/22/16.
//

#ifndef SIMDNA_DATASTRUCTURE_H
#define SIMDNA_DATASTRUCTURE_H

#include <stdio.h>
#include <iostream>
#include <unordered_map>
#include <vector>

#include "../lib/Vector3D.h"
#include "Configure.h"


/*
 * ID of each nucleotide residue {strand No., base No.}
 *
 *     default constructor assigns both strand and base number as -1
 *     therefore if there exists a base number as -1, meaning there is not such base or it hasn't defined
 * */
class ID {
    int _strandID;
    int _baseID;
public:

    ID() : _strandID{-1}, _baseID{-1} {}
    ID(int s, int b) : _strandID{s}, _baseID{b} {}

    int strandID() const { return _strandID; }
    int baseID() const { return _baseID; }

    bool operator<(const ID& rhs) const {
        return _strandID < rhs._strandID || (rhs._strandID >= _strandID && _baseID < rhs._baseID);
    }
    bool operator==(const ID& other) const {
        return (_strandID == other._strandID && _baseID == other._baseID);
    }

};

std::ostream& operator<<(std::ostream& os, const ID &idnum);

struct IDHasher {
    int operator()(const ID& k) const {
        return ((std::hash<int>()(k.strandID())
                 ^ (std::hash<int>()(k.baseID()) << 1)^(std::hash<int>()(k.baseID()) << 1)) >> 1);
    }
};

/*
 * Information (pair, coordinate, etc..) of a nucleotide
 *
 * */
class Nucleotide {
    ID _id;
    Vector3Dd _coordinate;
    ID _pairID;
    bool _isDS; // false if pair residue doesn't exist, _pairID == {-1, -1}
    bool _isBreak; // whether this residue belongs to helical break

public:
    Nucleotide() {}
    Nucleotide(ID IDnum, Vector3Dd coor, ID comID, bool isBreak) :
            _id{IDnum}, _pairID{comID}, _coordinate{coor},
            _isDS{ !((comID.strandID() == -1) && (comID.baseID() == -1))},
            _isBreak{ isBreak } {}


    ID id() const {return _id;}
    ID pairID() const {return _pairID;}
    Vector3Dd coordinate() const {return _coordinate;}
    bool withPair() const {return _isDS;}
    bool isBreak() const { return _isBreak;}

};

/*
 * Helical break pairs
 *
 *
 *
 * */
class HelicalBreakPair {

    bool _typeA;
    bool _typeB;
    bool _single;
    std::pair<ID, ID> _ids;


public:

    HelicalBreakPair() : _typeA(false), _typeB(false), _single(false) {}

    HelicalBreakPair(const std::pair<ID, ID> &_ids) :
            _ids(_ids), _typeA(false), _typeB(false), _single(false) {}

    HelicalBreakPair(const std::pair<ID, ID> &_ids, bool _typeA, bool _typeB, bool _single) :
            _typeA(_typeA), _typeB(_typeB), _single(_single), _ids(_ids) {}

    bool is_typeA() const { return _typeA; }
    void set_typeA(bool _typeA) { HelicalBreakPair::_typeA = _typeA; }
    bool is_typeB() const { return _typeB; }
    void set_typeB(bool _typeB) { HelicalBreakPair::_typeB = _typeB; }
    bool is_single() const { return _single; }
    void set_single(bool _single) { HelicalBreakPair::_single = _single; }
    const std::pair<ID, ID> &get_ids() const { return _ids; }

    static bool similar(const HelicalBreakPair& a, const HelicalBreakPair& b) {
        return a.is_typeA()==b.is_typeA()&&a.is_typeB()==b.is_typeB()&&
                a.is_single()==b.is_single();
    }
};




std::ostream& operator<<(std::ostream& os, const HelicalBreakPair & hbp);



class Node {
    int _num;
    int _type; // 3 types
    std::vector<std::pair<ID,ID>> _ids;
    double _mass;
    Vector3Dd _position;
    double _vdWradii;

public:


    Node() {}


    Node(int _type, const std::vector<std::pair<ID, ID>> &_ids,
         double _mass, const Vector3Dd &_position,
         double _vdWradii) : _type(_type), _ids(_ids), _mass(_mass),
                             _position(_position), _vdWradii(_vdWradii) {}

    void set_num(int _num) { Node::_num = _num; }
    int get_num() const { return _num; }
    int get_type() const { return _type; }
    void set_type(int _type) { Node::_type = _type; }
    const std::vector<std::pair<ID, ID>> &get_ids() const { return _ids; }
    void set_ids(const std::vector<std::pair<ID, ID>> &_ids) { Node::_ids = _ids; }
    double get_mass() const { return _mass; }
    void set_mass(double _mass) { Node::_mass = _mass; }
    const Vector3Dd &get_position() const { return _position; }
    void set_position(const Vector3Dd &_position) { Node::_position = _position; }
    double get_vdWradii() const { return _vdWradii; }
    void set_vdWradii(double _vdWradii) { Node::_vdWradii = _vdWradii; }


};

class Edge {
    std::pair<int, int> _endsNode;
    std::pair<int, int> _types;
    std::pair<ID, ID> _endsHB;

    bool _ds; // whether is double strand
    bool _crossover;
    std::vector<std::pair<ID,ID>> _ids;
    double _stretchConstant;

public:
    Edge() {}


    Edge(const std::pair<int, int> &_endsNode, const std::pair<int, int> &_types, const std::pair<ID, ID> &_endsHB,
         bool _crossover, const std::vector<std::pair<ID, ID>> &_ids)
            : _endsNode(_endsNode), _types(_types), _endsHB(_endsHB), _crossover(_crossover), _ids(_ids)
    {
        _ds = _types.first==1&&_types.second==1;
        _stretchConstant = _ds ? STRETCH_DS : STRETCH_SS;
    }

    Edge(const std::pair<ID, ID> &_endsHB, bool _crossover) :
            _endsHB(_endsHB), _crossover(_crossover) {}


    size_t length() const { return _ids.size(); }
    double get_stretchConstant() const { return _stretchConstant; }
    void set_stretchConstant(double _stretchConstant) { Edge::_stretchConstant = _stretchConstant; }
    const std::pair<int, int> &get_endsNode() const { return _endsNode; }
    void set_endsNode(const std::pair<int, int> &_endsNode) { Edge::_endsNode = _endsNode; }
    const std::pair<int, int> &get_types() const { return _types; }
    void set_types(const std::pair<int, int> &_types) { Edge::_types = _types; }
    const std::pair<ID, ID> &get_endsHB() const { return _endsHB; }
    void set_endsHB(const std::pair<ID, ID> &_endsHB) { Edge::_endsHB = _endsHB; }
    bool is_ds() const { return _ds; }
    void set_ds(bool _ds) { Edge::_ds = _ds; }
    bool is_crossover() const { return _crossover; }

    void set_crossover(bool _crossover) {
        Edge::_crossover = _crossover;
    }

    const std::vector<std::pair<ID, ID>> &get_ids() const {
        return _ids;
    }

    void set_ids(const std::vector<std::pair<ID, ID>> &_ids) {
        Edge::_ids = _ids;
    }

    void update() {
        ;

    }
};

/*
 * Mapping ID to abstract structure, i.e., helical break pairs, nodes, etc
 *      where Typename T is the abstracted structure
 *      all index starts from 1
 * */
template<typename T>
class Index {
    std::unordered_map<int, T> _member;
    std::unordered_map<ID, int, IDHasher> _index;
    int _size;

public:

    Index() : _size(0) {}
    int findIndexFromID(ID id) {return _index[id];}
    T findTypeFromIndex(int index) {return _member[index];}
    T findTypeFromID(ID id) { return _member[_index[id]];}
    int idExists(ID id) const {
        // if id exists, return the index No. associated to it, otherwise 0
        std::unordered_map<ID, int, IDHasher>::const_iterator a = _index.find(id);
        if ( a == _index.end()) return 0;
        else return a->second;
    }
    int size() {return _size;}
    std::unordered_map<int, T> member() const { return _member;}
    std::unordered_map<ID, int, IDHasher> index() const { return _index;}

    friend class HelicalBreakPairs;
    friend class Nodes;
    friend class Edges;


};


class HelicalBreakPairs : public Index<HelicalBreakPair> {

public:
    void insert(std::pair<ID, ID> basePair, bool a, bool b, bool c) {
        auto index = this->idExists(basePair.first);
        if (index) {
            if (a) _member[index].set_typeA(true);
            if (b) _member[index].set_typeB(true);
        }
        else {
            ++_size;
            _member[_size] = HelicalBreakPair{basePair, a, b, c};
            if (basePair.first.baseID()!=-1) _index[basePair.first] = _size;
            if (basePair.second.baseID()!=-1) _index[basePair.second] = _size;
        }
    }
    void print() const {
        std::cout << _size << std::endl;
        for (auto item : _member)
            std::cout << item.first << "\t" << item.second;
    }
};

class Nodes : public Index<Node> {
public:
    void insert(Node nd) {
        if (!this->idExists(nd.get_ids()[0].first)) {

            ++_size;
            if (_size==1099)
                std::cout << std::endl;
            nd.set_num(_size);
            _member[_size] = nd;
            for (const auto & item : nd.get_ids()) {
                if (item.first.baseID() != -1) _index[item.first] = _size;
                if (item.second.baseID() != -1) _index[item.second] = _size;
            }

        }
    }
};

class Edges : public Index<Edge> {
public:
    void insert(Edge eg) {
        ++_size;
        _member[_size] = eg;
        for (const auto & item : eg.get_ids()) {
            if (item.first.baseID() != -1) _index[item.first] = _size;
            if (item.second.baseID() != -1) _index[item.second] = _size;
        }

    }
    void update(Nodes &_nodes) {
        int a, b, c, d;
        ID id1, id2;
        Edge edge;
        for (auto && item : _member) {
            edge = item.second;
            a = _nodes.findIndexFromID(edge.get_endsHB().first);
            b = _nodes.findIndexFromID(edge.get_endsHB().second);
            edge.set_endsNode({a, b});
            c = _nodes.findTypeFromID(edge.get_endsHB().first).get_type();
            d = _nodes.findTypeFromID(edge.get_endsHB().second).get_type();
            edge.set_types({c, d});
        }
    }
};



class OrigamiGraph {
    Nodes _nodes;
    Edges _edges;
    std::vector<std::vector<int>> _origamiGraph; // adjacency list

public:
    void insertNode(Node nd) { _nodes.insert(nd); }
    void insertEdge(Edge eg, int c1, int c2);
    bool findEdge(int id1, int id2) const;
    int findNodeType(ID id);
    int findNodeType(int a);
    const std::vector<std::vector<int>>& showGraph() const { return _origamiGraph; }

    int findNodeNum(ID id);
    int howManyNodes() { return _nodes.size();}
    Vector3Dd nodeCenter(int a) { return _nodes.findTypeFromIndex(a).get_position();}
    void printAllNodes() {
        for (auto && item :_nodes.member())
            for (auto && item1 : item.second.get_ids())
                std::cout << item1.first << "\t" << item1.second << std::endl;
    }
    void resizeGraph() {
        _origamiGraph.resize(_nodes.size()+1);
    }
    void printID(int a, int b, int item1) {
        ID id;
        id = _nodes.findTypeFromIndex(a).get_ids()[item1].first;
//        if (id.baseID() == 2533) {
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].second << std::endl;
////            std::cout << _nodes.findTypeFromIndex(a).get_ids().size() << std::endl;
//            std::cout << "+++++++++++++++++++" << std::endl;
//        }
        std::cout << id << std::endl;
//        if (id.baseID() == 2501) std::cout << a << "\t" << b<< "1"<< std::endl;

        id = _nodes.findTypeFromIndex(a).get_ids()[item1].second;
//        if (id.baseID() == 2533) {
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].second << std::endl;
////            std::cout << _nodes.findTypeFromIndex(a).get_ids().size() << std::endl;
//            std::cout << "+++++++++++++++++++" << std::endl;
//        }
        std::cout << id << std::endl;
//        if (id.baseID() == 2501) std::cout << a << "\t" << b<< "2"<< std::endl;

        id = _nodes.findTypeFromIndex(b).get_ids()[item1].first;
        std::cout << id << std::endl;
//        if (id.baseID() == 2533) {
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].second << std::endl;
////            std::cout << a << "\t" << b << std::endl;
//
////            std::cout << _nodes.findTypeFromIndex(a).get_ids().size() << std::endl;
//
//            std::cout << "+++++++++++++++++++" << std::endl;
//        }

        id = _nodes.findTypeFromIndex(b).get_ids()[item1].second;
//        if (id.baseID() == 2533) {
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[1].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(a).get_ids()[0].second << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].first << std::endl;
//            std::cout << _nodes.findTypeFromIndex(b).get_ids()[0].second << std::endl;
////            std::cout << _nodes.findTypeFromIndex(a).get_ids().size() << std::endl;
//            std::cout << "+++++++++++++++++++" << std::endl;
//        }
        std::cout << id << std::endl;
//        if (id.baseID() == 2501) std::cout << a << "\t" << b<< "4"<< std::endl;


    }
    void howManyFourWays() {
        Edge edge;
        int sum = 0;
        int a, b;
        ID id;
        for (const auto& item : _edges.member()) {
            edge = item.second;
            if (edge.get_types().first == 1 && edge.get_types().second == 1 && edge.is_crossover()) {
                ++sum;

                a = edge.get_endsNode().first;
                b = edge.get_endsNode().second;

                for (int item1=0; item1 < _nodes.findTypeFromIndex(a).get_ids().size(); ++item1) {
//                    printID(a, b, item1);
                }

            }
        }
        std::cout << sum << std::endl;

    }

};


#endif //SIMDNA_DATASTRUCTURE_H

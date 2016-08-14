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

    bool operator<(const ID &rhs) const {
        return _strandID < rhs._strandID || (rhs._strandID >= _strandID && _baseID < rhs._baseID);
    }
    bool operator==(const ID &other) const {
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

struct pairhash {
public:
    template <typename T, typename U>
    std::size_t operator()(const std::pair<T, U> &x) const
    {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }
};

/*
 * Nucleotide information (ID, pair residue ID, coordinate, etc..)
 *
 * */
class Nucleotide {
    ID _id;
    Vector3Dd _coordinate;
    ID _pairID;
    bool _isDS; // false if pair residue doesn't exist, _pairID == {-1, -1}
    bool _isBreak; // whether this residue belongs to axial discontinuity 

public:
    Nucleotide() {}
    Nucleotide(const ID &IDnum, const Vector3Dd &coor, const ID &comID, bool isBreak) :
            _id{IDnum}, _pairID{comID}, _coordinate{coor},
            _isDS{ !((comID.strandID() == -1) && (comID.baseID() == -1))},
            _isBreak{ isBreak } {}

    const ID &id() const {return _id;}
    const ID &pairID() const {return _pairID;}
    const Vector3Dd &coordinate() const {return _coordinate;}
    bool withPair() const {return _isDS;}
    bool isDiscont() const { return _isBreak;}

};

/*
 * Axial discontinuity
 * A base pair belong to axial discontinuity
 * Type A: Neighbor to another discontinuity not by a crossover
 * Type B: Neighbor to another discontinuity by a crossover
 * Single: doesn't have pair
 *
 * For example: (assume all "=" duplex, "-" single strand)
 * ==i=c======a=b======d
 *            | |      |
 * k-j=e======f=g======h
 * a, b, f, and g are both type A and B; d and h are type B; i, j and c are type A; k is single
 * */
class AxialDiscontinuity {

    bool _typeA;
    bool _typeB;
    bool _single;
    std::pair<ID, ID> _ids;

public:

    AxialDiscontinuity() : _typeA(false), _typeB(false), _single(false) {}
    AxialDiscontinuity(const std::pair<ID, ID> &_ids) :
            _ids(_ids), _typeA(false), _typeB(false), _single(false) {}
    AxialDiscontinuity(const std::pair<ID, ID> &_ids, bool _typeA, bool _typeB, bool _single) :
            _typeA(_typeA), _typeB(_typeB), _single(_single), _ids(_ids) {}


    bool is_typeA() const { return _typeA; }
    void set_typeA(bool _typeA) { AxialDiscontinuity::_typeA = _typeA; }
    bool is_typeB() const { return _typeB; }
    void set_typeB(bool _typeB) { AxialDiscontinuity::_typeB = _typeB; }
    bool is_single() const { return _single; }
    const std::pair<ID, ID> &get_ids() const { return _ids; }

    static bool similar(const AxialDiscontinuity &a, const AxialDiscontinuity &b) {
        return a.is_typeA()==b.is_typeA()&&a.is_typeB()==b.is_typeB()&&
                a.is_single()==b.is_single();
    }
};

std::ostream& operator<<(std::ostream& os, const AxialDiscontinuity & hbp);


/*
 * Basic unit used in the simulation
 *
 * */
class Node {
    int _num;
    int _type; // 3 types
    std::vector<std::pair<ID,ID>> _ids;
    double _mass;
    Vector3Dd _position;
    double _vdWradii;
    bool _inHJ; // whether in a Holliday junction

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
    bool is_inHJ() const { return _inHJ; }
    void set_inHJ(bool _inHJ) { Node::_inHJ = _inHJ; }
};



/*
 * Basic unit used in the simulation
 *
 *
 * */
class Edge {
    std::pair<int, int> _endsNode;
    std::pair<int, int> _types;
    std::pair<ID, ID> _endsHB;

    bool _ds; // whether is double strand
    bool _crossover; // whether is a crossover
    bool _isHJ; // whether is a Holliday junction
    bool _isStackJunc; // whether is in stacked junctions
    std::vector<std::pair<ID,ID>> _ids;
    double _stretchConstant;

public:
    Edge() {}
    Edge(const std::pair<int, int> &_endsNode, const std::pair<int, int> &_types,
         const std::pair<ID, ID> &_endsHB, bool _crossover,
         const std::vector<std::pair<ID, ID>> &_ids)
            : _endsNode(_endsNode), _types(_types), _endsHB(_endsHB), _crossover(_crossover),
              _ids(_ids), _isHJ(false), _isStackJunc(false)
    {
        _ds = _crossover ? _types.first == 1 && _types.second == 1 : _types.first != 3 && _types.second != 3;
        _stretchConstant = _ds ? STRETCH_DS : STRETCH_SS;
    }
    Edge(const std::pair<int, int> &_endsNode) // initialize stack only
            : _endsNode(_endsNode), _types({2,2}), _endsHB({ID{}, ID{}}), _ds(true), _crossover(false),
               _isHJ(false), _isStackJunc(true), _ids({}), _stretchConstant(STRETCH_SS/100) {}

    size_t length() const { return _ids.size(); }
    double get_stretchConstant() const { return _stretchConstant; }
    const std::pair<int, int> &get_endsNode() const { return _endsNode; }
    void set_endsNode(const std::pair<int, int> &_endsNode) { Edge::_endsNode = _endsNode; }
    const std::pair<int, int> &get_types() const { return _types; }
    void set_types(const std::pair<int, int> &_types) { Edge::_types = _types; }
    const std::pair<ID, ID> &get_endsHB() const { return _endsHB; }
    bool is_ds() const { return _ds; }
    bool is_crossover() const { return _crossover; }
    const std::vector<std::pair<ID, ID>> &get_ids() const { return _ids; }
    bool is_isHJ() const { return _isHJ; }
    void set_isHJ(bool _isHJ) { Edge::_isHJ = _isHJ; }
    bool is_isStackJunc() const { return _isStackJunc; }
};

/**
 * Maps ID with an arbitrary class, i.e., Edge, Node, etc.
 *      where Typename T is the type of the arbitrary class,
 *      all index starts from 1.
 */
template<typename T>
class Index {
    std::unordered_map<int, T> _member;
    std::unordered_map<ID, int, IDHasher> _index;
    int _size;

public:

    Index() : _size(0) {}

    int &findIndexFromID(ID id) {return _index[id];}
    T &findTypeFromIndex(int index) {return _member[index];}
    T &findTypeFromID(ID id) { return _member[_index[id]];}

    int size() {return _size;}
    const std::unordered_map<int, T> &member() const { return _member;}
    std::unordered_map<int, T> &updateMember() { return _member;}
    const std::unordered_map<ID, int, IDHasher> &index() const { return _index;}

    int idExists(ID id) const {
        // if id exists, return the index No., otherwise 0
        std::unordered_map<ID, int, IDHasher>::const_iterator a = _index.find(id);
        if ( a == _index.end()) return 0;
        else return a->second;
    }

    friend class AxialDiscons;
    friend class Nodes;
    friend class Edges;

};


class AxialDiscons : public Index<AxialDiscontinuity> {

public:
    void insert(std::pair<ID, ID> basePair, bool a, bool b, bool c) {
        auto index = this->idExists(basePair.first);
        if (index) {
            if (a) _member[index].set_typeA(true);
            if (b) _member[index].set_typeB(true);
        }
        else {
            ++_size;
            _member[_size] = AxialDiscontinuity{basePair, a, b, c};
            if (basePair.first.baseID()!=-1) _index[basePair.first] = _size;
            if (basePair.second.baseID()!=-1) _index[basePair.second] = _size;
        }
    }
};

class Nodes : public Index<Node> {
public:
    void insert(Node nd);
    int insertVir(Node nd);
};

class Edges : public Index<Edge> {
    std::unordered_map<std::pair<int, int>, int, pairhash> _indexEdge;
public:
    void insert(Edge eg) {
        ++_size;
        _member[_size] = eg;
        _indexEdge[eg.get_endsNode()] = _size;
        for (const auto &item : eg.get_ids()) {
            if (item.first.baseID() != -1) _index[item.first] = _size;
            if (item.second.baseID() != -1) _index[item.second] = _size;
        }
    }
    const Edge &findEdgeFromEnds(int a, int b) const {


        auto v = _indexEdge.find({a, b});
        if ( v == _indexEdge.end()) return _member.at(_indexEdge.at({b, a}));
        else return _member.at(v->second);

    }

};



class OrigamiGraph {
    Nodes _nodes;
    Edges _edges;
    std::vector<std::vector<int>> _origamiGraph;
    int _numHJ; // how many Holliday junctions
    std::vector<Edge> _HJ;
    std::vector<Edge> _crossovers;

public:

    void resizeGraph() { _origamiGraph.resize(_nodes.size()+1); }
    void insertNode(Node nd) { _nodes.insert(nd); }
    int insertVir(Node nd) { return _nodes.insertVir(nd);  }
    void insertEdge(Edge eg, int c1, int c2);
    bool findEdge(int id1, int id2) const;
    int findNodeType(ID id);
    int findNodeType(int a);
    int findNodeNum(ID id);
    const std::vector<std::pair<ID, ID>> &findIDsfromIndex(int a);
    const std::vector<std::vector<int>>& showGraph() const { return _origamiGraph; }
    int howManyNodes() { return _nodes.size();}
    Vector3Dd nodeCenter(int a) { return _nodes.findTypeFromIndex(a).get_position();}
    const Edges &get_edges() const { return _edges; }
    const Nodes &get_nodes() const { return _nodes; }
    const std::vector<Edge> &get_HJ() const { return _HJ; }
    const std::vector<Edge> &get_crossovers() const { return _crossovers; }
    const std::vector<int> connectFrom(int a) const { return _origamiGraph[a]; }
    void repCon(int a, int b, int c);

    int howMany4Ways() const { return _numHJ; }
    const Edge &findEdgeFromEnds(int a, int b) const {
        return _edges.findEdgeFromEnds(a, b);
    }
    const Node &findNodeFromNum(int a) { return _nodes.findTypeFromIndex(a);}

    void printAllNodes() {
        for (auto && item :_nodes.member())
            for (auto && item1 : item.second.get_ids())
                std::cout << item1.first << "\t" << item1.second << std::endl;
    }

/*
 * Find all Holliday junctions, complete _nodes and _edges
 * */
    void findHollidayJ();

};

class VirtualSite {
    int virNum;
    Vector3Dd coordinate;
    int end1;
    int end2;
    double weight1;
    double weight2;
    
public:
    VirtualSite(int virNum, const Vector3Dd &coordinate, int end1, int end2, double weight1, double weight2) : virNum(
            virNum), coordinate(coordinate), end1(end1), end2(end2), weight1(weight1), weight2(weight2) { }
    
    double getWeight2() const { return weight2; }
    int getVirNum() const { return virNum; }
    const Vector3Dd &getCoordinate() const { return coordinate; }
    int getEnd1() const { return end1; }
    int getEnd2() const { return end2; }
    double getWeight1() const { return weight1; }
};
#endif //SIMDNA_DATASTRUCTURE_H

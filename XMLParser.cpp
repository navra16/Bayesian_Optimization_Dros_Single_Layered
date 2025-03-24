#include "XMLParser.h"
#include <iostream>
#include <cstdio>

bool XMLParser::parseFile(const std::string& filename, SystemBuilder& builder) {
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_file(filename.c_str());
    
    std::cout<<"data struct file name = "<< filename.c_str()<<std::endl;
    if (!result) {
        std::cerr << "XMLParser: Error parsing file: " << result.description() << std::endl;
        return false;
    }
    
    pugi::xml_node root = doc.child("data");
    if (!root) {
        std::cerr << "XMLParser: No <data> root found." << std::endl;
        return false;
    }
    
    // ---------------------------
    // Process <settings> node
    pugi::xml_node settings = root.child("settings");
    if (settings) {
        if (auto p = settings.child("Tau")) {
            builder.defaultTau = p.text().as_double();
            std::cout << "Setting tau: " << builder.defaultTau << std::endl;
        }
        if (auto p = settings.child("KBT")) {
            builder.defaultKBT = p.text().as_double();
            std::cout << "Setting kbt: " << builder.defaultKBT << std::endl;
        }
        if (auto p = settings.child("Linear_Const")) {
            builder.defaultLinear_Const = p.text().as_double();
            std::cout << "Setting linear const: " << builder.defaultLinear_Const << std::endl;
        }
        if (auto p = settings.child("Area_Const")) {
            builder.defaultArea_Const = p.text().as_double();
            std::cout << "Setting area const: " << builder.defaultArea_Const << std::endl;
        }
        if (auto p = settings.child("Bend_Const")) {
            builder.defaultBending_Const = p.text().as_double();
            std::cout << "Setting bending const: " << builder.defaultBending_Const << std::endl;
        }
        if (auto p = settings.child("LJ_Eps")) {
            builder.defaultLJ_Eps = p.text().as_double();
            std::cout << "Setting lj eps: " << builder.defaultLJ_Eps << std::endl;
        }
        if (auto p = settings.child("LJ_Rmin")) {
            builder.defaultLJ_Rmin = p.text().as_double();
            std::cout << "Setting lj rmin: " << builder.defaultLJ_Rmin << std::endl;
        }
        if (auto p = settings.child("LJ_Rmax")) {
            builder.defaultLJ_Rmax = p.text().as_double();
            std::cout << "Setting lj rmax: " << builder.defaultLJ_Rmax << std::endl;
        }
        if (auto p = settings.child("LJ_Const")) {
            builder.defaultLJ_Const = p.text().as_double();
            std::cout << "Setting lj const: " << builder.defaultLJ_Const << std::endl;
        }
        if (auto p = settings.child("LJ_X")) {
            builder.defaultLJ_X = p.text().as_double();
            std::cout << "Setting lj x: " << builder.defaultLJ_X << std::endl;
        }
        if (auto p = settings.child("LJ_Y")) {
            builder.defaultLJ_Y = p.text().as_double();
            std::cout << "Setting lj y: " << builder.defaultLJ_Y << std::endl;
        }
        if (auto p = settings.child("LJ_Z")) {
            builder.defaultLJ_Z = p.text().as_double();
            std::cout << "Setting lj z: " << builder.defaultLJ_Z << std::endl;
        }
        // (If additional settings are added in the future, they can be parsed here.)
    }
    
    // ---------------------------
    // Process <nodes> section
    pugi::xml_node nodes = root.child("nodes");
    if (nodes) {
        for (pugi::xml_node node = nodes.child("node"); node; node = node.next_sibling("node")) {
            double x, y, z;
            const char* text = node.text().as_string();
            if (3 != sscanf(text, "%lf %lf %lf", &x, &y, &z)) {
                std::cerr << "XMLParser: parse node error" << std::endl;
                return false;
            }
            builder.addNode(x, y, z);
        }
        std::cout << "Nodes parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Process <edgeinfos> section
    pugi::xml_node edgeinfos = root.child("edgeinfos");
    if (edgeinfos) {
        for (pugi::xml_node edge = edgeinfos.child("edgeinfo"); edge; edge = edge.next_sibling("edgeinfo")) {
            unsigned int from, to;
            double restLength;
            int count = sscanf(edge.text().as_string(), "%u %u %lf", &from, &to, &restLength);
            if (count == 3) {
                builder.addEdge(from - 1, to - 1, restLength);
            } else if (count == 2) {
                builder.addEdge(from - 1, to - 1, builder.defaultEdgeEq);
            } else {
                std::cerr << "XMLParser: parse edge error" << std::endl;
                return false;
            }
        }
        std::cout << "Edges parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Process <elems> section (elements/triangles)
    pugi::xml_node elems = root.child("elems");
    if (elems) {
        for (pugi::xml_node elem = elems.child("elem"); elem; elem = elem.next_sibling("elem")) {
            unsigned int a, b, c;
            if (3 != sscanf(elem.text().as_string(), "%u %u %u", &a, &b, &c)) {
                std::cerr << "XMLParser: parse elem error" << std::endl;
                return false;
            }
            builder.addElement(a - 1, b - 1, c - 1);
        }
        std::cout << "Elements parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Process <elem2edges> mapping, if present.
    pugi::xml_node elem2edges = root.child("elem2edges");
    if (elem2edges) {
        for (pugi::xml_node mapping = elem2edges.child("elem2edge"); mapping; mapping = mapping.next_sibling("elem2edge")) {
            unsigned int e1, e2, e3;
            if (3 != sscanf(mapping.text().as_string(), "%u %u %u", &e1, &e2, &e3)) {
                std::cerr << "XMLParser: parse elem2edge error" << std::endl;
                return false;
            }
            builder.addElement2Edge(e1 - 1, e2 - 1, e3 - 1);
        }
        std::cout << "Element-to-edge mappings parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Process <edge2elems> mapping, if present.
    pugi::xml_node edge2elems = root.child("edge2elems");
    if (edge2elems) {
        for (pugi::xml_node mapping = edge2elems.child("edge2elem"); mapping; mapping = mapping.next_sibling("edge2elem")) {
            unsigned int elemID, edgeID;
            if (2 != sscanf(mapping.text().as_string(), "%u %u", &elemID, &edgeID)) {
                std::cerr << "XMLParser: parse edge2elem error" << std::endl;
                return false;
            }
            builder.addEdge2Elem(elemID - 1, edgeID - 1);
        }
        std::cout << "Edge-to-element mappings parsed successfully." << std::endl;
    }
    
    // ---------------------------
    // Optionally process other sections (capsidnodes, fixed nodes, nndatas, etc.)
    pugi::xml_node capsidnodes = root.child("capsidnodes");
    if (capsidnodes) {
        for (pugi::xml_node node = capsidnodes.child("node"); node; node = node.next_sibling("node")) {
            double x, y, z;
            if (3 != sscanf(node.text().as_string(), "%lf %lf %lf", &x, &y, &z)) {
                std::cerr << "XMLParser: parse capsid node error" << std::endl;
                return false;
            }
            builder.addCapsidNode(x, y, z);
        }
        std::cout << "Capsid nodes parsed successfully." << std::endl;
    }
    
    pugi::xml_node fixnodes = root.child("fixed");
    if (fixnodes) {
        for (pugi::xml_node node = fixnodes.child("node"); node; node = node.next_sibling("node")) {
            int id = node.text().as_int();
            builder.fixNodes(id - 1);
        }
        std::cout << "Fixed nodes parsed successfully." << std::endl;
    }
    
    pugi::xml_node nndatas = root.child("nndatas");
    if (nndatas) {
        for (pugi::xml_node nndata = nndatas.child("nndata"); nndata; nndata = nndata.next_sibling("nndata")) {
            double d1, d2, d3, d4, d5, d6, d7, d8, d9;
            if (9 != sscanf(nndata.text().as_string(), "%lf %lf %lf %lf %lf %lf %lf %lf %lf",
                            &d1, &d2, &d3, &d4, &d5, &d6, &d7, &d8, &d9)) {
                std::cerr << "XMLParser: parse nndata error" << std::endl;
                return false;
            }
            builder.addNndata(d1, d2, d3, d4, d5, d6, d7, d8, d9);
        }
        std::cout << "nndata parsed successfully." << std::endl;
    }
    
    return true;
}

#include "BarnesHut.h"
#include <cassert>

class Node{
            public:
                long double mass = 0; // total mass of node
                vector<long double> com; // centre of mass of node
                vector<GravitationalBody*> particles; //ptrs to particles inside the node
                vector<long double> centre; // centre of node
                long double s; // side length
                vector<Node*> child; // ptrs to child nodes
                Node(long double x, long double y, long double z, long double side){
                    centre = {x,y,z};
                    s = side;
                }
};

class BarnesHutMethods {
    private:
        void Divide(Node* node){
            long double x = node->centre[0];
            long double y = node->centre[1];
            long double z = node->centre[2];
            long double s = node->s;
            Node* btr = new Node(x+s/4,y+s/4,z-s/4,s/2); // b-bottom t-top l-left r-right
            Node* btl = new Node(x-s/4,y+s/4,z-s/4,s/2);
            Node* bbl = new Node(x-s/4,y-s/4,z-s/4,s/2);
            Node* bbr = new Node(x+s/4,y-s/4,z-s/4,s/2);
            Node* ttr = new Node(x+s/4,y+s/4,z+s/4,s/2);
            Node* ttl = new Node(x-s/4,y+s/4,z+s/4,s/2);
            Node* tbl = new Node(x-s/4,y-s/4,z+s/4,s/2);
            Node* tbr = new Node(x+s/4,y-s/4,z+s/4,s/2);
            node->child = {btr,btl,bbl,bbr,ttr,ttl,tbl,tbr};
        }

        int whichOctant(Node* node, GravitationalBody& body){
            long double xr = body.position[0] - node->centre[0];
            long double yr = body.position[1] - node->centre[1];
            long double zr = body.position[2] - node->centre[2];
            if (xr>0 && yr>0 && zr<0){return 0;}
            else if (xr<0 && yr>0 && zr<0){return 1;}
            else if (xr<0 && yr<0 && zr<0){return 2;}
            else if (xr>0 && yr<0 && zr<0){return 3;}
            else if (xr>0 && yr>0 && zr>0){return 4;}
            else if (xr<0 && yr>0 && zr>0){return 5;}
            else if (xr<0 && yr<0 && zr>0){return 6;}
            else {return 7;}
        }

        void Octinsert(Node* node, GravitationalBody& body){
            if (node->particles.size()>1){
                node->com[0] = (node->mass*node->com[0] + body.mass*body.position[0])/(node->mass+body.mass);
                node->com[1] = (node->mass*node->com[1] + body.mass*body.position[1])/(node->mass+body.mass);
                node->com[2] = (node->mass*node->com[2] + body.mass*body.position[2])/(node->mass+body.mass);
                node->mass += body.mass;
                node->particles.push_back(&body);
                Octinsert(node->child[whichOctant(node,body)],body);
            }
            else if (node->particles.size()==1){
                Divide(node);
                Octinsert(node->child[whichOctant(node,*(node->particles[0]))],*(node->particles[0]));
                node->com[0] = (node->mass*node->com[0] + body.mass*body.position[0])/(node->mass+body.mass);
                node->com[1] = (node->mass*node->com[1] + body.mass*body.position[1])/(node->mass+body.mass);
                node->com[2] = (node->mass*node->com[2] + body.mass*body.position[2])/(node->mass+body.mass);
                node->mass += body.mass;
                node->particles.push_back(&body);
                Octinsert(node->child[whichOctant(node,body)],body);
            }
            else {
                node->com = body.position;
                node->mass = body.mass;
                node->particles = {&body};
            }
        }

        void DeleteUnusedNodes(Node* root){
            for (auto it=root->child.begin();it!=root->child.end();){
                if ((*it)->particles.size()==0){
                    delete *it;
                    it = root->child.erase(it);
                }
                else if ((*it)->particles.size()==1){
                    it++;
                }
                else {
                    DeleteUnusedNodes(*it);
                    it++;
                }
            }
        }

        void OcttreeBuild(Node* root, vector<GravitationalBody>& o){
            for (int i=0;i<o.size();i++){
                Octinsert(root,o[i]);
            }
            DeleteUnusedNodes(root);
        }

        void DeleteOcttree(Node* root){
            for (auto it=root->child.begin();it!=root->child.end();it++){
                if ((*it)->particles.size()==1){
                    delete *it;
                }
                else {
                    DeleteOcttree(*it);
                }
            }
            delete root;
        }

        vector<long double> TreeForce(Node* node, GravitationalBody& body){
            long double fx=0, fy=0, fz=0;
            if (node->particles.size()==1){
                if (node->particles[0]!=&body){
                    long double r2 = pow((node->com[0]-body.position[0]),2)+pow((node->com[1]-body.position[1]),2)+pow((node->com[2]-body.position[2]),2);
                    fx = G*body.mass*node->mass*(node->com[0]-body.position[0])/pow(Q_rsqrt(r2),3);
                    fy = G*body.mass*node->mass*(node->com[1]-body.position[1])/pow(Q_rsqrt(r2),3);
                    fz = G*body.mass*node->mass*(node->com[2]-body.position[2])/pow(Q_rsqrt(r2),3);                    
                }
            }
            else {
                long double r = Q_rsqrt(pow((node->com[0]-body.position[0]),2)+pow((node->com[1]-body.position[1]),2)+pow((node->com[2]-body.position[2]),2));
                long double theta = 0.5;
                if (node->s/r <theta){
                    fx = G*body.mass*node->mass*(node->com[0]-body.position[0])/pow(r,3);
                    fy = G*body.mass*node->mass*(node->com[1]-body.position[1])/pow(r,3);
                    fz = G*body.mass*node->mass*(node->com[2]-body.position[2])/pow(r,3);
                }
                else {
                    for (auto it=node->child.begin();it!=node->child.end();it++){
                        vector<long double> ChildForce = TreeForce(*it,body);
                        fx+=ChildForce[0];
                        fy+=ChildForce[1];
                        fz+=ChildForce[2];
                    }
                }
            }
            return {fx,fy,fz};
        }

    public:
        vector<long double> RootCenter;
        long double RootSide;
        BarnesHutMethods(vector<long double> center, long double side){
            RootCenter = center;
            RootSide = side;
        }
        vector<vector<long double>> BarnesHutForces(vector<GravitationalBody>& o){
            Node* root = new Node(RootCenter[0],RootCenter[1],RootCenter[2],RootSide);
            OcttreeBuild(root,o);
            vector<long double> fx(o.size(),0), fy(o.size(),0), fz(o.size(),0);
            for (int i=0;i<o.size();i++){
                vector<long double> ForceOnBody = TreeForce(root,o[i]);
                fx[i] = ForceOnBody[0];
                fy[i] = ForceOnBody[1];
                fz[i] = ForceOnBody[2];
            }
            DeleteOcttree(root);
            return {fx,fy,fz};
        }
};

BarnesHut::BarnesHut(GravitationalSystem& s):ForceCalculator(s){
    BarnesHutMethods TreeObj({0, 0, 0}, 100);                         //need to change based on length scale of the problem(rename later)
    forces = TreeObj.BarnesHutForces(s.bodies);
}

valtype BarnesHut::getForce(const int i, const int coordType){
    assert(forces.size() == 3);
    assert(forces[i].size() == forces[0].size());                                  //Segmentiation errors, Check
    //cout << i << endl;
    return forces[coordType][i];
}

#include <iostream>
#include <queue>
#include <set>
#include <unordered_map>
#include <vector>

#include "../geometry/halfedge.h"
#include "debug.h"

/* Note on local operation return types:

    The local operations all return a std::optional<T> type. This is used so that your
    implementation can signify that it does not want to perform the operation for
    whatever reason (e.g. you don't want to allow the user to erase the last vertex).

    An optional can have two values: std::nullopt, or a value of the type it is
    parameterized on. In this way, it's similar to a pointer, but has two advantages:
    the value it holds need not be allocated elsewhere, and it provides an API that
    forces the user to check if it is null before using the value.

    In your implementation, if you have successfully performed the operation, you can
    simply return the required reference:

            ... collapse the edge ...
            return collapsed_vertex_ref;

    And if you wish to deny the operation, you can return the null optional:

            return std::nullopt;

    Note that the stubs below all reject their duties by returning the null optional.
*/

/*
    This method should replace the given vertex and all its neighboring
    edges and faces with a single face, returning the new face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_vertex(Halfedge_Mesh::VertexRef v) {
    std::cout << __FUNCTION__ << std::endl;
    auto h = v->halfedge();
    auto face = this->new_face();
    auto h_boundary = h;
    do {
        auto h_twin = h->twin();
        auto h_1 = h->next();
        auto h_v = h_1->vertex();

        // The edge and half edge will be deleted
        this->erase(h);
        this->erase(h_twin);
        this->erase(h->edge());
        this->erase(h_twin->face());

        // If the edge deleted, make sure the boundary vertex has right outgoing half edge
        if(h_v->halfedge() == h_twin) h_v->halfedge() = h_1;

        h_boundary = h_twin->next();
        do {
            h_boundary = h_boundary->next();
            // Set all boundary to the new face
            h_boundary->face() = face;
        } while(h_boundary->next() != h_twin);

        // The boundary should connect to boundary half edge, for the inner edge has been deleted.
        h_boundary->next() = h_1;

        h = h_twin->next();
    } while(h != v->halfedge());
    face->halfedge() = h_boundary; // Assign a half boundary to new face

    this->erase(v);
    return face;
}

/*
    This method should erase the given edge and return an iterator to the
    merged face.
 */
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::erase_edge(Halfedge_Mesh::EdgeRef e) {
    std::cout << __FUNCTION__ << std::endl;

    // Handle the extreme case when flag ==1
    if(e->halfedge()->next() == e->halfedge()->twin()) e->halfedge() = e->halfedge()->twin();

    auto h = e->halfedge();
    auto h_twin = h->twin();

    int flag = 0;
    if(h_twin->next() == h) flag = 1;

    auto v = h->vertex();
    auto v_bair = h_twin->vertex();

    auto face = this->new_face();
    // 1. register all halfedges on the old face

    if(flag == 1) {
        // This is a extreme case, when the h and h_twin have the same face.
        do {
            h = h->next();
            // register the face to the half edge.
            h->face() = face;
        } while(h->next() != h_twin);

        // Make a loop in new face half edges
        h->next() = e->halfedge()->next();

        this->erase(v);
        v_bair->halfedge() = h->next();
    } else {
        // 1. register all halfedges on the old face
        do {
            h = h->next();
            // register the face to the half edge.
            h->face() = face;
        } while(h->next() != e->halfedge());

        // Make a loop in new face half edges
        h->next() = h_twin->next();

        // Set the old vertex outgoing edge.
        v->halfedge() = h->next();
        // 2. register all halfedges on another side old face
        do {
            h_twin = h_twin->next();
            h_twin->face() = face;
        } while(h_twin->next() != e->halfedge()->twin());
        h_twin->next() = e->halfedge()->next(); // close the loop
        v_bair->halfedge() = h_twin->next();
        // Set the old vertex outgoing edge.
        v_bair->halfedge() = h_twin->next();
    }

    // register a halfedge to the newface
    face->halfedge() = h->next();

    // delete all components related to this edge
    this->erase(e->halfedge());
    this->erase(e->halfedge()->twin());
    this->erase(e->halfedge()->face());
    this->erase(e->halfedge()->twin()->face());
    this->erase(e);

    return face;
}

/*
    This method should collapse the given edge and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_edge(Halfedge_Mesh::EdgeRef e) {
    std::cout << __FUNCTION__ << std::endl;

    // This will deal with a manifold with arbitrary edges
    // For collapsing we don't need to add any new elements, just delete some and reassign some
    // edges and faces
    auto h_db = e->halfedge();
    auto h_ba = h_db->next();
    auto h_ax = h_ba->next();
    auto h_xd = h_ax;
    while(h_xd->next() != h_db) {
        h_xd = h_xd->next();
    }
    bool bf1_triangle = false;
    if(h_xd == h_ax) bf1_triangle = true;

    auto h_bd = h_db->twin();
    auto h_dw = h_bd->next();
    auto h_wc = h_dw;
    while(h_wc->next()->next() != h_bd) {
        h_wc = h_wc->next();
    }
    auto h_cb = h_wc->next();
    bool bf2_triangle = false;
    if(h_wc == h_dw) bf2_triangle = true;

    auto h_ab = h_ba->twin();
    auto h_be = h_ab->next();
    auto h_ey = h_be->next();
    auto h_ya = h_ey;
    while(h_ya->next() != h_ab) {
        h_ya = h_ya->next();
    }

    auto h_bc = h_cb->twin();
    auto h_cz = h_bc->next();
    auto h_zf = h_cz;
    while(h_zf->next()->next() != h_bc) {
        h_zf = h_zf->next();
    }
    auto h_fb = h_zf->next();

    auto v_a = h_ax->vertex();
    auto v_b = h_bd->vertex();
    auto v_c = h_cz->vertex();
    auto v_d = h_db->vertex();

    auto f_1 = h_db->face();
    auto f_2 = h_bd->face();
    auto f_3 = h_be->face();
    auto f_4 = h_fb->face();

    // We certainly need to erase collapse edge and corresponding vertex and half edges
    this->erase(e);
    this->erase(v_b);
    this->erase(h_db);
    this->erase(h_bd);

    // Reassign affected vertices, face and half edges.
    v_d->pos = (v_d->pos + v_b->pos) / 2;
    if(v_d->halfedge() == h_db) v_d->halfedge() = h_dw;

    if(bf1_triangle) {
        // f_1 is a triangle, we need to eliminate the face and edge ab to avoid duplicate
        this->erase(h_ba->edge());
        this->erase(f_1);
        this->erase(h_ba);
        this->erase(h_ab);

        if(v_a->halfedge() == h_ab) v_a->halfedge() = h_ax;
        if(f_3->halfedge() == h_ab) f_3->halfedge() = h_be;

        h_ya->set_neighbors(h_ax, h_ya->twin(), h_ya->vertex(), h_ya->edge(), f_3);
        h_xd->set_neighbors(h_be, h_xd->twin(), v_a, h_xd->edge(), f_3);
    } else {
        if(f_1->halfedge() == h_db) f_1->halfedge() = h_ba;
        h_ba->vertex() = v_d;
        h_xd->next() = h_ba;
    }

    if(bf2_triangle) {
        // f_2 is a triangle, we need to elimate the face and edge bc to avoid duplicate
        this->erase(h_bc->edge());
        this->erase(h_bc);
        this->erase(h_cb);
        this->erase(f_2);

        if(v_c->halfedge() == h_cb) v_c->halfedge() = h_cz;
        if(f_4->halfedge() == h_bc) f_4->halfedge() = h_fb;

        h_dw->set_neighbors(h_cz, h_dw->twin(), v_d, h_dw->edge(), f_4);
        h_fb->set_neighbors(h_dw, h_fb->twin(), h_fb->vertex(), h_fb->edge(), f_4);
    } else {
        if(f_2->halfedge() == h_bd) f_2->halfedge() = h_cb;
        h_bc->vertex() = v_d;
        h_cb->next() = h_dw;
    }

    // All half edge outgoing from vertex b should assignt to d now
    auto h = h_be;
    do {
        h->vertex() = v_d;
        h = h->twin()->next();
    } while(h != h_fb->twin());
    h->vertex() = v_d;

    return v_d;
}

/*
    This method should collapse the given face and return an iterator to
    the new vertex created by the collapse.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::collapse_face(Halfedge_Mesh::FaceRef f) {
    std::cout << __FUNCTION__ << std::endl;
    (void)f;
    return std::nullopt;
}

/*
    This method should flip the given edge and return an iterator to the
    flipped edge.
*/
std::optional<Halfedge_Mesh::EdgeRef> Halfedge_Mesh::flip_edge(Halfedge_Mesh::EdgeRef e) {
    std::cout << __FUNCTION__ << std::endl;

    // Get all half edges associate with flip operations
    auto h_bd = e->halfedge();
    auto h_da = h_bd->next();
    auto h_ax = h_da->next();
    auto h_xb = h_ax;
    while(h_xb->next() != h_bd) {
        h_xb = h_xb->next();
    }

    auto h_db = h_bd->twin();
    auto h_bc = h_db->next();
    auto h_cy = h_bc->next();
    auto h_yd = h_cy;
    while(h_yd->next() != h_db) {
        h_yd = h_yd->next();
    }

    // The two vertices connected edge e
    auto v_a = h_ax->vertex();
    auto v_b = h_bc->vertex();
    auto v_c = h_cy->vertex();
    auto v_d = h_da->vertex();

    // The two face on this edge e
    auto f_1 = h_bd->face();
    auto f_2 = h_db->face();

    // The two vertices going to connected new edge

    // reset the connected vertex with proper outgoing halfedge
    if(v_b->halfedge() == h_bd) v_b->halfedge() = h_bc;
    if(v_d->halfedge() == h_db) v_d->halfedge() = h_da;

    //    v_a->halfedge() = h;
    //    v_c->halfedge() = h_twin;

    // reset the information on the face
    if(f_1->halfedge() == h_da) f_1->halfedge() = h_bd;
    if(f_2->halfedge() == h_bc) f_2->halfedge() = h_db;

    // reset the information on all associate halfedges
    h_bd->set_neighbors(h_ax, h_db, v_c, e, f_1);
    h_xb->set_neighbors(h_bc, h_xb->twin(), h_xb->vertex(), h_xb->edge(), f_1);
    h_bc->set_neighbors(h_bd, h_bc->twin(), v_b, h_bc->edge(), f_1);

    h_db->set_neighbors(h_cy, h_bd, v_a, e, f_2);
    h_yd->set_neighbors(h_da, h_yd->twin(), h_yd->vertex(), h_yd->edge(), f_2);
    h_da->set_neighbors(h_db, h_da->twin(), v_d, h_da->edge(), f_2);

    return e;
}

/*
    This method should split the given edge and return an iterator to the
    newly inserted vertex. The halfedge of this vertex should point along
    the edge that was split, rather than the new edges.
*/
std::optional<Halfedge_Mesh::VertexRef> Halfedge_Mesh::split_edge(Halfedge_Mesh::EdgeRef e) {
    std::cout << __FUNCTION__ << std::endl;

    // Get the associate vertices and half edges
    auto h_ac = e->halfedge();
    auto h_cd = h_ac->next();
    auto h_da = h_cd->next();

    auto h_ca = h_ac->twin();
    auto h_ab = h_ca->next();
    auto h_bc = h_ab->next();

    // vertex a,c are the vertex associate to the edge,
    auto v_a = h_ac->vertex();
    auto v_c = h_ca->vertex();

    // for the edge is going to be deleted, if anybody associate with this edge, reassign them
    if(v_a->halfedge() == h_ac) v_a->halfedge() = h_ab;
    if(v_c->halfedge() == h_ca) v_c->halfedge() = h_cd;

    auto v_d = h_da->vertex();
    auto v_b = h_bc->vertex();

    // The half edges and faces associate this edge shall be deleted
    this->erase(h_ac);
    this->erase(h_ca);
    this->erase(e);
    this->erase(h_ac->face());
    this->erase(h_ca->face());

    // Split the edge will create 4 new faces, 4 edges, 8 halfedges and 1 vertex
    auto h_am = this->new_halfedge();
    auto h_ma = this->new_halfedge();
    auto h_bm = this->new_halfedge();
    auto h_mb = this->new_halfedge();
    auto h_cm = this->new_halfedge();
    auto h_mc = this->new_halfedge();
    auto h_dm = this->new_halfedge();
    auto h_md = this->new_halfedge();

    // Create a new vertex and associate with new half edge
    auto v_m = this->new_vertex();
    v_m->halfedge() = h_mc;
    v_m->pos = (v_a->pos + v_c->pos) / 2;

    // Create 4 new faces and asscoiate with proper half edges
    auto f_1 = this->new_face(); // associate with triangle dam
    f_1->halfedge() = h_da;
    auto f_2 = this->new_face(); // associate with triangle abm
    f_2->halfedge() = h_ab;
    auto f_3 = this->new_face(); // associate with triangle bcm
    f_3->halfedge() = h_bc;
    auto f_4 = this->new_face(); // associate with triangle cdm
    f_4->halfedge() = h_cd;

    // Create 4 new edges and associate with proper half edges
    auto e_1 = this->new_edge(); // associate with edge am
    e_1->halfedge() = h_am;
    auto e_2 = this->new_edge(); // associate with edge bm
    e_2->halfedge() = h_mb;
    auto e_3 = this->new_edge(); // associate with edge cm
    e_3->halfedge() = h_mc;
    auto e_4 = this->new_edge(); // associate with edge dm
    e_4->halfedge() = h_dm;

    // Associate proper info to each half edge
    h_am->set_neighbors(h_md, h_ma, v_a, e_1, f_1);
    h_ma->set_neighbors(h_ab, h_am, v_m, e_1, f_2);
    h_bm->set_neighbors(h_ma, h_mb, v_b, e_2, f_2);
    h_mb->set_neighbors(h_bc, h_bm, v_m, e_2, f_3);
    h_cm->set_neighbors(h_mb, h_mc, v_c, e_3, f_3);
    h_mc->set_neighbors(h_cd, h_cm, v_m, e_3, f_4);
    h_dm->set_neighbors(h_mc, h_md, v_d, e_4, f_4);
    h_md->set_neighbors(h_da, h_dm, v_m, e_4, f_1);

    // Loop the triangles on the old half edges
    h_da->next() = h_am;
    h_da->face() = f_1;
    h_ab->next() = h_bm;
    h_ab->face() = f_2;
    h_bc->next() = h_cm;
    h_bc->face() = f_3;
    h_cd->next() = h_dm;
    h_cd->face() = f_4;

    return v_m;
}

/* Note on the beveling process:

    Each of the bevel_vertex, bevel_edge, and bevel_face functions do not represent
    a full bevel operation. Instead, they should update the _connectivity_ of
    the mesh, _not_ the positions of newly created vertices. In fact, you should set
    the positions of new vertices to be exactly the same as wherever they "started from."

    When you click on a mesh element while in bevel mode, one of those three functions
    is called. But, because you may then adjust the distance/offset of the newly
    beveled face, we need another method of updating the positions of the new vertices.

    This is where bevel_vertex_positions, bevel_edge_positions, and
    bevel_face_positions come in: these functions are called repeatedly as you
    move your mouse, the position of which determins the normal and tangent offset
    parameters. These functions are also passed an array of the original vertex
    positions: for  bevel_vertex, it has one element, the original vertex position,
    for bevel_edge,  two for the two vertices, and for bevel_face, it has the original
    position of each vertex in halfedge order. You should use these positions, as well
    as the normal and tangent offset fields to assign positions to the new vertices.

    Finally, note that the normal and tangent offsets are not relative values - you
    should compute a particular new position from them, not a delta to apply.
*/

/*
    This method should replace the vertex v with a face, corresponding to
    a bevel operation. It should return the new face.  NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_vertex_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_vertex(Halfedge_Mesh::VertexRef v) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."
    std::cout << __FUNCTION__ << std::endl;

    auto h = v->halfedge();
    auto f_new = this->new_face();

    std::vector<Halfedge_Mesh::VertexRef> vs;
    std::vector<Halfedge_Mesh::HalfedgeRef> hs, h_twins;
    std::vector<Halfedge_Mesh::EdgeRef> es;
    // just prelocate do not do anything operation yet
    do {
        auto h_twin = h->twin();
        // Associate all new elements
        vs.push_back(this->new_vertex());
        hs.push_back(this->new_halfedge());
        h_twins.push_back(this->new_halfedge());
        es.push_back(this->new_edge());

        h = h_twin->next();
    } while(h != v->halfedge());

    int i = 0;
    for(i = 0, h = v->halfedge(); i < vs.size(); ++i) {
        auto h_new = hs[i];
        auto h_twin_new = h_twins[i];
        auto h_twin_new_p = h_twins[(i - 1 + vs.size()) % vs.size()];

        auto v_0 = vs[i];
        auto v_1 = vs[(i + 1) % vs.size()];
        auto e_new = es[i];

        auto h_twin = h->twin();
        auto h_end = h;
        auto f = h->face();
        while(h_end->next() != h) {
            h_end = h_end->next();
        }

        v_0->halfedge() = h_new;
        v_0->pos = v->pos;
        e_new->halfedge() = h_new;
        h_end->next() = h_new;
        h_new->set_neighbors(h, h_twin_new, v_0, e_new, f);
        h->vertex() = v_1;
        h_twin_new->set_neighbors(h_twin_new_p, h_new, v_1, e_new, f_new);

        h = h_twin->next();
    }
    this->erase(v);
    f_new->halfedge() = h_twins[0];

    return f_new;
}

/*
    This method should replace the edge e with a face, corresponding to a
    bevel operation. It should return the new face. NOTE: This method is
    responsible for updating the *connectivity* of the mesh only---it does not
    need to update the vertex positions.  These positions will be updated in
    Halfedge_Mesh::bevel_edge_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_edge(Halfedge_Mesh::EdgeRef e) {
    std::cout << __FUNCTION__ << std::endl;

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    auto h_1 = e->halfedge();
    auto h_2 = h_1->twin();
    auto f_1 = h_1->face();
    auto f_2 = h_2->face();
    auto ve_1 = h_1->vertex();
    auto ve_2 = h_2->vertex();
    auto f_new = this->new_face();

    this->erase(e);
    this->erase(h_1);
    this->erase(h_2);
    this->erase(ve_1);
    this->erase(ve_2);

    std::vector<Halfedge_Mesh::VertexRef> vs;
    std::vector<Halfedge_Mesh::HalfedgeRef> hs, h_twins;
    std::vector<Halfedge_Mesh::EdgeRef> es;
    // Prelocate some elements, start from h_1 face side, until reach h_2
    h_1 = h_1->next();
    do {
        auto h_twin = h_1->twin();
        // Associate all new elements
        vs.push_back(this->new_vertex());
        hs.push_back(this->new_halfedge());
        h_twins.push_back(this->new_halfedge());
        es.push_back(this->new_edge());

        h_1 = h_twin->next();
    } while(h_1 != h_2);
    h_1 = e->halfedge(); // reset h_1

    // now start from h_2 face side until reach h_1;

    h_2 = h_2->next();
    do {
        auto h_twin = h_2->twin();
        // Associate all new elements
        vs.push_back(this->new_vertex());
        hs.push_back(this->new_halfedge());
        h_twins.push_back(this->new_halfedge());
        es.push_back(this->new_edge());

        h_2 = h_twin->next();
    } while(h_2 != h_1);
    h_2 = e->halfedge()->twin(); // reset h_2

    int i = 0;
    auto h = h_1;
    for(i = 0, h = e->halfedge()->next(); i < vs.size(); ++i) {
        auto h_new = hs[i];
        auto h_twin_new = h_twins[i];
        auto h_twin_new_p = h_twins[(i - 1 + vs.size()) % vs.size()];

        auto v_0 = vs[i];
        auto v_1 = vs[(i + 1) % vs.size()];
        auto e_new = es[i];

        if(h == h_2) h = h_2->next();

        auto h_twin = h->twin();
        auto h_end = h;
        auto f = h->face();
        if(h->face() == f_1 || h->face() == f_2) {
            h->face()->halfedge() = h;
            while(h_end->next()->next() != h) {
                h_end = h_end->next();
            }
        } else {
            while(h_end->next() != h) {
                h_end = h_end->next();
            }
        }

        v_0->halfedge() = h_new;
        if(h->face() == f_1) {
            v_0->pos = h_1->vertex()->pos;
        } else if(h->face() == f_2) {
            v_0->pos = h_2->vertex()->pos;
        } else {
            v_0->pos = h->vertex()->pos;
        }
        e_new->halfedge() = h_new;
        h_end->next() = h_new;
        h_new->set_neighbors(h, h_twin_new, v_0, e_new, f);
        h->vertex() = v_1;
        h_twin_new->set_neighbors(h_twin_new_p, h_new, v_1, e_new, f_new);

        h = h_twin->next();
    }

    f_new->halfedge() = h_twins[0];

    return f_new;
}

/*
    This method should replace the face f with an additional, inset face
    (and ring of faces around it), corresponding to a bevel operation. It
    should return the new face.  NOTE: This method is responsible for updating
    the *connectivity* of the mesh only---it does not need to update the vertex
    positions. These positions will be updated in
    Halfedge_Mesh::bevel_face_positions (which you also have to
    implement!)
*/
std::optional<Halfedge_Mesh::FaceRef> Halfedge_Mesh::bevel_face(Halfedge_Mesh::FaceRef f) {

    // Reminder: You should set the positions of new vertices (v->pos) to be exactly
    // the same as wherever they "started from."

    (void)f;
    return std::nullopt;
}

/*
    Compute new vertex positions for the vertices of the beveled vertex.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the original vertex position and its associated outgoing edge
    to compute a new vertex position along the outgoing edge.
*/
void Halfedge_Mesh::bevel_vertex_positions(const std::vector<Vec3>& start_positions,
                                           Halfedge_Mesh::FaceRef face, float tangent_offset) {
    std::cout << __FUNCTION__ << std::endl;

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    for(int i = 0; i < new_halfedges.size(); ++i) {
        h = new_halfedges[i];
        auto h_twin = h->twin();

        auto v_target = h->vertex();
        auto v_n = h_twin->next()->next()->vertex();

        Vec3 dir = (v_n->pos - start_positions[i]).normalize();
        v_target->pos = start_positions[i] + tangent_offset * dir;
    }
}

/*
    Compute new vertex positions for the vertices of the beveled edge.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the orig array) to compute an offset vertex position.

    Note that there is a 1-to-1 correspondence between halfedges in
    newHalfedges and vertex positions
    in orig.  So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vector3D pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_edge_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset) {
    std::cout << __FUNCTION__ << std::endl;

    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    for(int i = 0; i < new_halfedges.size(); ++i) {
        h = new_halfedges[i];
        auto h_twin = h->twin();

        auto v_target = h->vertex();
        auto v_n = h_twin->next()->next()->vertex();

        Vec3 dir = (v_n->pos - start_positions[i]).normalize();
        v_target->pos = start_positions[i] + tangent_offset * dir;
    }
}

/*
    Compute new vertex positions for the vertices of the beveled face.

    These vertices can be accessed via new_halfedges[i]->vertex()->pos for
    i = 1, ..., new_halfedges.size()-1.

    The basic strategy here is to loop over the list of outgoing halfedges,
    and use the preceding and next vertex position from the original mesh
    (in the start_positions array) to compute an offset vertex
    position.

    Note that there is a 1-to-1 correspondence between halfedges in
    new_halfedges and vertex positions
    in orig. So, you can write loops of the form

    for(size_t i = 0; i < new_halfedges.size(); i++)
    {
            Vec3 pi = start_positions[i]; // get the original vertex
            position corresponding to vertex i
    }
*/
void Halfedge_Mesh::bevel_face_positions(const std::vector<Vec3>& start_positions,
                                         Halfedge_Mesh::FaceRef face, float tangent_offset,
                                         float normal_offset) {

    if(flip_orientation) normal_offset = -normal_offset;
    std::vector<HalfedgeRef> new_halfedges;
    auto h = face->halfedge();
    do {
        new_halfedges.push_back(h);
        h = h->next();
    } while(h != face->halfedge());

    (void)new_halfedges;
    (void)start_positions;
    (void)face;
    (void)tangent_offset;
    (void)normal_offset;
}

/*
    Splits all non-triangular faces into triangles.
*/
void Halfedge_Mesh::triangulate() {

    // For each face...
}

/* Note on the quad subdivision process:

        Unlike the local mesh operations (like bevel or edge flip), we will perform
        subdivision by splitting *all* faces into quads "simultaneously."  Rather
        than operating directly on the halfedge data structure (which as you've
        seen is quite difficult to maintain!) we are going to do something a bit nicer:
           1. Create a raw list of vertex positions and faces (rather than a full-
              blown halfedge mesh).
           2. Build a new halfedge mesh from these lists, replacing the old one.
        Sometimes rebuilding a data structure from scratch is simpler (and even
        more efficient) than incrementally modifying the existing one.  These steps are
        detailed below.

  Step I: Compute the vertex positions for the subdivided mesh.
        Here we're going to do something a little bit strange: since we will
        have one vertex in the subdivided mesh for each vertex, edge, and face in
        the original mesh, we can nicely store the new vertex *positions* as
        attributes on vertices, edges, and faces of the original mesh. These positions
        can then be conveniently copied into the new, subdivided mesh.
        This is what you will implement in linear_subdivide_positions() and
        catmullclark_subdivide_positions().

  Steps II-IV are provided (see Halfedge_Mesh::subdivide()), but are still detailed
  here:

  Step II: Assign a unique index (starting at 0) to each vertex, edge, and
        face in the original mesh. These indices will be the indices of the
        vertices in the new (subdivided mesh).  They do not have to be assigned
        in any particular order, so long as no index is shared by more than one
        mesh element, and the total number of indices is equal to V+E+F, i.e.,
        the total number of vertices plus edges plus faces in the original mesh.
        Basically we just need a one-to-one mapping between original mesh elements
        and subdivided mesh vertices.

  Step III: Build a list of quads in the new (subdivided) mesh, as tuples of
        the element indices defined above. In other words, each new quad should be
        of the form (i,j,k,l), where i,j,k and l are four of the indices stored on
        our original mesh elements.  Note that it is essential to get the orientation
        right here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces
        should circulate in the same direction as old faces (think about the right-hand
        rule).

  Step IV: Pass the list of vertices and quads to a routine that clears
        the internal data for this halfedge mesh, and builds new halfedge data from
        scratch, using the two lists.
*/

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    simple linear interpolation, e.g., the edge midpoints and face
    centroids.
*/
void Halfedge_Mesh::linear_subdivide_positions() {

    // For each vertex, assign Vertex::new_pos to
    // its original position, Vertex::pos.

    // For each edge, assign the midpoint of the two original
    // positions to Edge::new_pos.

    // For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::new_pos. Note
    // that in general, NOT all faces will be triangles!
}

/*
    Compute new vertex positions for a mesh that splits each polygon
    into quads (by inserting a vertex at the face midpoint and each
    of the edge midpoints).  The new vertex positions will be stored
    in the members Vertex::new_pos, Edge::new_pos, and
    Face::new_pos.  The values of the positions are based on
    the Catmull-Clark rules for subdivision.

    Note: this will only be called on meshes without boundary
*/
void Halfedge_Mesh::catmullclark_subdivide_positions() {

    // The implementation for this routine should be
    // a lot like Halfedge_Mesh:linear_subdivide_positions:(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    // Faces

    // Edges

    // Vertices
}

/*
        This routine should increase the number of triangles in the mesh
        using Loop subdivision. Note: this is will only be called on triangle meshes.
*/
void Halfedge_Mesh::loop_subdivide() {

    // Compute new positions for all the vertices in the input mesh, using
    // the Loop subdivision rule, and store them in Vertex::new_pos.
    // -> At this point, we also want to mark each vertex as being a vertex of the
    //    original mesh. Use Vertex::is_new for this.
    // -> Next, compute the updated vertex positions associated with edges, and
    //    store it in Edge::new_pos.
    // -> Next, we're going to split every edge in the mesh, in any order.  For
    //    future reference, we're also going to store some information about which
    //    subdivided edges come from splitting an edge in the original mesh, and
    //    which edges are new, by setting the flat Edge::is_new. Note that in this
    //    loop, we only want to iterate over edges of the original mesh.
    //    Otherwise, we'll end up splitting edges that we just split (and the
    //    loop will never end!)
    // -> Now flip any new edge that connects an old and new vertex.
    // -> Finally, copy the new vertex positions into final Vertex::pos.

    // Each vertex and edge of the original surface can be associated with a
    // vertex in the new (subdivided) surface.
    // Therefore, our strategy for computing the subdivided vertex locations is to
    // *first* compute the new positions
    // using the connectivity of the original (coarse) mesh; navigating this mesh
    // will be much easier than navigating
    // the new subdivided (fine) mesh, which has more elements to traverse.  We
    // will then assign vertex positions in
    // the new mesh based on the values we computed for the original mesh.

    // Compute updated positions for all the vertices in the original mesh, using
    // the Loop subdivision rule.

    // Next, compute the updated vertex positions associated with edges.

    // Next, we're going to split every edge in the mesh, in any order. For
    // future reference, we're also going to store some information about which
    // subdivided edges come from splitting an edge in the original mesh, and
    // which edges are new.
    // In this loop, we only want to iterate over edges of the original
    // mesh---otherwise, we'll end up splitting edges that we just split (and
    // the loop will never end!)

    // Finally, flip any new edge that connects an old and new vertex.

    // Copy the updated vertex positions to the subdivided mesh.
}

/*
    Isotropic remeshing. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if this is not a triangle mesh)
*/
bool Halfedge_Mesh::isotropic_remesh() {

    // Compute the mean edge length.
    // Repeat the four main steps for 5 or 6 iterations
    // -> Split edges much longer than the target length (being careful about
    //    how the loop is written!)
    // -> Collapse edges much shorter than the target length.  Here we need to
    //    be EXTRA careful about advancing the loop, because many edges may have
    //    been destroyed by a collapse (which ones?)
    // -> Now flip each edge if it improves vertex degree
    // -> Finally, apply some tangential smoothing to the vertex positions

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

/* Helper type for quadric simplification */
struct Edge_Record {
    Edge_Record() {
    }
    Edge_Record(std::unordered_map<Halfedge_Mesh::VertexRef, Mat4>& vertex_quadrics,
                Halfedge_Mesh::EdgeRef e)
        : edge(e) {

        // Compute the combined quadric from the edge endpoints.
        // -> Build the 3x3 linear system whose solution minimizes the quadric error
        //    associated with these two endpoints.
        // -> Use this system to solve for the optimal position, and store it in
        //    Edge_Record::optimal.
        // -> Also store the cost associated with collapsing this edge in
        //    Edge_Record::cost.
    }
    Halfedge_Mesh::EdgeRef edge;
    Vec3 optimal;
    float cost;
};

/* Comparison operator for Edge_Records so std::set will properly order them */
bool operator<(const Edge_Record& r1, const Edge_Record& r2) {
    if(r1.cost != r2.cost) {
        return r1.cost < r2.cost;
    }
    Halfedge_Mesh::EdgeRef e1 = r1.edge;
    Halfedge_Mesh::EdgeRef e2 = r2.edge;
    return &*e1 < &*e2;
}

/** Helper type for quadric simplification
 *
 * A PQueue is a minimum-priority queue that
 * allows elements to be both inserted and removed from the
 * queue.  Together, one can easily change the priority of
 * an item by removing it, and re-inserting the same item
 * but with a different priority.  A priority queue, for
 * those who don't remember or haven't seen it before, is a
 * data structure that always keeps track of the item with
 * the smallest priority or "score," even as new elements
 * are inserted and removed.  Priority queues are often an
 * essential component of greedy algorithms, where one wants
 * to iteratively operate on the current "best" element.
 *
 * PQueue is templated on the type T of the object
 * being queued.  For this reason, T must define a comparison
 * operator of the form
 *
 *    bool operator<( const T& t1, const T& t2 )
 *
 * which returns true if and only if t1 is considered to have a
 * lower priority than t2.
 *
 * Basic use of a PQueue might look
 * something like this:
 *
 *    // initialize an empty queue
 *    PQueue<myItemType> queue;
 *
 *    // add some items (which we assume have been created
 *    // elsewhere, each of which has its priority stored as
 *    // some kind of internal member variable)
 *    queue.insert( item1 );
 *    queue.insert( item2 );
 *    queue.insert( item3 );
 *
 *    // get the highest priority item currently in the queue
 *    myItemType highestPriorityItem = queue.top();
 *
 *    // remove the highest priority item, automatically
 *    // promoting the next-highest priority item to the top
 *    queue.pop();
 *
 *    myItemType nextHighestPriorityItem = queue.top();
 *
 *    // Etc.
 *
 *    // We can also remove an item, making sure it is no
 *    // longer in the queue (note that this item may already
 *    // have been removed, if it was the 1st or 2nd-highest
 *    // priority item!)
 *    queue.remove( item2 );
 *
 */
template<class T> struct PQueue {
    void insert(const T& item) {
        queue.insert(item);
    }
    void remove(const T& item) {
        if(queue.find(item) != queue.end()) {
            queue.erase(item);
        }
    }
    const T& top(void) const {
        return *(queue.begin());
    }
    void pop(void) {
        queue.erase(queue.begin());
    }
    size_t size() {
        return queue.size();
    }

    std::set<T> queue;
};

/*
    Mesh simplification. Note that this function returns success in a similar
    manner to the local operations, except with only a boolean value.
    (e.g. you may want to return false if you can't simplify the mesh any
    further without destroying it.)
*/
bool Halfedge_Mesh::simplify() {

    std::unordered_map<VertexRef, Mat4> vertex_quadrics;
    std::unordered_map<FaceRef, Mat4> face_quadrics;
    std::unordered_map<EdgeRef, Edge_Record> edge_records;
    PQueue<Edge_Record> edge_queue;

    // Compute initial quadrics for each face by simply writing the plane equation
    // for the face in homogeneous coordinates. These quadrics should be stored
    // in face_quadrics
    // -> Compute an initial quadric for each vertex as the sum of the quadrics
    //    associated with the incident faces, storing it in vertex_quadrics
    // -> Build a priority queue of edges according to their quadric error cost,
    //    i.e., by building an Edge_Record for each edge and sticking it in the
    //    queue. You may want to use the above PQueue<Edge_Record> for this.
    // -> Until we reach the target edge budget, collapse the best edge. Remember
    //    to remove from the queue any edge that touches the collapsing edge
    //    BEFORE it gets collapsed, and add back into the queue any edge touching
    //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
    //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
    //    top of the queue.

    // Note: if you erase elements in a local operation, they will not be actually deleted
    // until do_erase or validate are called. This is to facilitate checking
    // for dangling references to elements that will be erased.
    // The rest of the codebase will automatically call validate() after each op,
    // but here simply calling collapse_edge() will not erase the elements.
    // You should use collapse_edge_erase() instead for the desired behavior.

    return false;
}

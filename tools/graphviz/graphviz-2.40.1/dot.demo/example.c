#include <gvc.h>

#define NO_LAYOUT_OR_RENDERING

int main(int argc, char **argv)
{
    Agraph_t *g;
    Agnode_t *n, *m;
    Agedge_t *e;

#ifndef NO_LAYOUT_OR_RENDERING
    /* set up a graphviz context - but only once even for multiple graphs */
    static GVC_t *gvc;

    if (!gvc)
    	gvc = gvContext();
#endif

    /* Create a simple digraph */
    g = agopen("g", AGDIGRAPH);
    n = agnode(g, "n");
    m = agnode(g, "m");
    e = agedge(g, n, m);

    /* Set an attribute - in this case one that affects the visible rendering */
    agsafeset(n, "color", "red", "");
    
#ifdef NO_LAYOUT_OR_RENDERING
    /* Just write the graph without layout */
    agwrite(g, stdout);
#else
    /* Use the directed graph layout engine */
    gvLayout(gvc, g, "dot");

    /* Output in .dot format */
    gvRender(gvc, g, "dot", stdout);

    gvFreeLayout(gvc, g);
#endif

    agclose(g);

    return 0;
}

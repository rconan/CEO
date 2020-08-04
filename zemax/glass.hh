
typedef struct _GlassEntry {
    char*	comment;
    char*	name;
    int		formula;
    double	MIL;
    double	Nd;
    double	Vd;
    int		exclude;
    int		status;
    double	TCE;
    double	TCE100300;
    double	density;
    double	dPgF;
    int		ignthermal;
    double	c[10];
    double	D0, D1, D2, E0, E1, Ltk, temp;
    double	ymin, ymax;
} GlassEntry;

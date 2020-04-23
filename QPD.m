(* ::Package:: *)

BeginPackage[ "QPD`"]

AutoRotate::usage = "AutoRotate[gr_Graphics3D, rate] takes a Graphics3D object and creates a dynamic rotating view that completes one rotation in about *rate* seconds."

CrystalPlot::usage = "Make a plot for a crystal."

ImportCIF::usage="ImportCIF[file] imports a crystal structure from the given CIF file, and returns a dataset that includes the space group of the structure, the name of the file used in querying the function, information about the constituent atoms, a mathematica-module object as parsed directly by the standard Import function. The returned crystal structure and four versions of the possible interpretations that the CIF file can have: the asymmetric ( and the complete cell, and the chemical and non-chemical types for each. These four versions may not actually be different depending on the input file. The final structure in the list is always the complete cell with chemical atom types

returning multiple versions of the structure as explained below. \
(ImportCIF can only handle files with one structure, but files with multiple structures can easily be split with any text editor.)

CIF files do not store all atoms in a structure explicitly. Instead, they contain only those atoms that are not connected to one another by \
symmetry operations, or the atoms of the asymmetric unit cell. The file then also contains all symmetry operations of the cell, so that the \
remaining atoms can be constructed. This is what ImportCIF will do. Additionally, the atom types stored in a CIF file may not be actual chemical \
types, they can also be a concatenation of the chemical type and an integer unique to the atom. Once the other atoms have been added by symmetry, \
the types will be unique to the atom's symmetry position. Sometimes these \"non-chemical\" types are helpful to have, although in most cases you \
will want the chemical types. To allow access to all of this data, ImportCIF returns four structure datasets for each crystal structure: the \
asymmetric and the complete cell, and the chemical and non-chemical types for each. (These four datasets may or may not actually be different \
from one another depending on the input file, but you always get all four of them.) The comment in the last entry of each structure will let you \
know which version of the structure you're dealing with. The final structure in the list is always the complete cell with chemical atom types, \
so Last[ImportCIF[file]] will give you the most physically meaningful version of the structure.";

Begin["Private`"]

numericKeysAll = {"_atom_site_fract_x", "_atom_site_fract_y", 
    "_atom_site_fract_z",
    "_cell_length_a", "_cell_length_b", "_cell_length_c", 
    "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma", 
    "_cell_volume"};
chosenProperties = {EntityProperty["Element", "AtomicRadius"], 
    EntityProperty["Element", "CovalentRadius"], 
    EntityProperty["Element", "ElectronConfiguration"], 
    EntityProperty["Element", "ElectronConfigurationString"], 
    EntityProperty["Element", "IconColor"], 
    EntityProperty["Element", "MostCommonOxidationStates"], 
    EntityProperty["Element", "Name"], 
    EntityProperty["Element", "RefractiveIndex"], 
    EntityProperty["Element", "SpaceGroup"], 
    EntityProperty["Element", "StableIsotopes"], 
    EntityProperty["Element", "VanDerWaalsRadius"]};

CrystalPlot[lattvec_, coord_, conf_, label_, sysdim_, addQ_, atomcol_, radii_, linestyle_]:=

    Module[{position, coordP, confP, new, tuples},   (*Module[{position, coordP, confP, new, tuples, bonds}, *)
    (*periodic repetitions: *)
    coordP=Flatten[Table[(# + {a, b, c}) &/@ coord, {a, 0, sysdim[[1]]-1}, 
                                                  {b, 0, sysdim[[2]]-1},
                                                  {c, 0, sysdim[[3]]-1}], 3];
    confP=Flatten[ConstantArray[conf, sysdim[[1]] * sysdim[[2]] * sysdim[[3]]]];

    (*add periodic duplicates if desired: *)
    If[TrueQ[addQ], 
        new={coordP, confP}\[Transpose];
        (*along the x-axis*)
        new=Join[new, 
            (# + {{sysdim[[1]], 0, 0}, 0}) &/@ Select[new, Abs[#[[1, 1]]]<0.01&], 
            (# + {{-sysdim[[1]], 0, 0}, 0}) &/@ Select[new, Abs[#[[1, 1]]]>0.99*sysdim[[1]]&]];
        (*along the y-axis*)
        new=Join[new, 
            (# + {{0, sysdim[[2]], 0}, 0}) &/@ Select[new, Abs[#[[1, 2]]]<0.01&], 
            (# + {{0, -sysdim[[2]], 0}, 0}) &/@ Select[new, Abs[#[[1, 2]]]>0.99*sysdim[[2]]&]];
        (*along the z-axis*)
        new=Join[new, 
            (# + {{0, 0, sysdim[[3]]}, 0}) &/@ Select[new, Abs[#[[1, 3]]]<0.01&], 
            (# + {{0, 0, -sysdim[[3]]}, 0}) &/@ Select[new, Abs[#[[1, 3]]]>0.99*sysdim[[3]]&]];
        coordP=new\[Transpose][[1]];
        confP=new\[Transpose][[2]];
    ];

    (*
    (*all pairs of atoms to be connected by bonds: *)
    tuples=Select[Tuples[Range[Length[confP]], 2],
                  (TrueQ[#[[1]] <# [[2]]] 
                   && Norm[(coordP[[#]].lattvec)[[1]] - (coordP[[#]].lattvec)[[2]]] < bondstyle[[1]]) &];
    bonds=coordP[[#]].lattvec&/@tuples;*)

    Graphics3D[{
    Specularity[White, 100], 
    (*cell outlines: *)
    linestyle[[2]], 
    Tube[#.lattvec, linestyle[[3]]]&/@Which[
        (*unit cell style: *)
        TrueQ[linestyle[[1]]==1], 
        {{{0, 0, #}, {1, 0, #}, {1, 1, #}, {0, 1, #}, {0, 0, #}} &/@ {0, 1}, 
        {{0, 0, #} &/@ {0, 1}, 
        {1, 0, #} &/@ {0, 1}, 
        {1, 1, #} &/@ {0, 1}, 
        {0, 1, #} &/@ {0, 1}}}, 
        (*periodic repetition style: *)
        TrueQ[linestyle[[1]]==2], 
        Join[Flatten[Table[{{0, #, z}, {sysdim[[1]], #, z}} &/@ Range[0, sysdim[[2]]], {z, 0, sysdim[[3]]}], 1], 
             Flatten[Table[{{#, y, 0}, {#, y, sysdim[[3]]}} &/@ Range[0, sysdim[[1]]], {y, 0, sysdim[[2]]}], 1], 
             Flatten[Table[{{x, 0, #}, {x, sysdim[[2]], #}} &/@ Range[0, sysdim[[3]]], {x, 0, sysdim[[1]]}], 1]], 
        (*supercell style: *)
        TrueQ[linestyle[[1]]==3], 
        #.(sysdim*IdentityMatrix[3]) &/@ Join[
            {{0, 0, #}, {1, 0, #}, {1, 1, #}, {0, 1, #}, {0, 0, #}} &/@ {0, 1}, 
            {{0, 0, #} &/@ {0, 1}, 
            {1, 0, #} &/@ {0, 1}, 
            {1, 1, #} &/@ {0, 1}, 
            {0, 1, #} &/@ {0, 1}}], 
        (*none: *)
        True, {}], 

    (*arrows for lattice vectors, unit cell/long/none: *)
    (*linestyle[[2]],*)
    Which[
        MemberQ[{1, 2}, linestyle[[1]]], 
        {Arrowheads[.07], MapThread[{#2,Arrow[Tube[{{0, 0, 0}, #1}, 3*linestyle[[3]]]]}&,{lattvec, {Blue, Red, Green}}]}, 
        TrueQ[linestyle[[1]]==3], 
        {Arrowheads[.07], Arrow[Tube[{{0, 0, 0}, #}, linestyle[[3]]]]&/@(lattvec*sysdim)}, 
        True, {}
        ], 

    (*atoms: *)
    {atomcol[[confP[[#]]]], HatchShading[], Sphere[coordP[[#]].lattvec, radii[[confP[[#]]]]]} &/@Range[Length[confP]]

    (*
    (*bonds: *)
    Which[
        (*bonds in atom colours: *)
        TrueQ[bondstyle[[2]]==1], 
        Table[{HatchShading[0.8, atomcol[[confP[[tuples[[#, ii]]]]]]], 
              Tube[{bonds[[#, ii]], Total[bonds[[#]]]*.5}, 
                bondstyle[[4]]]}&/@Range[Length[bonds]], 
              {ii, 1, 2}], 
        (*bonds in vertex colours interpolated between atom colours: *)
        TrueQ[bondstyle[[2]]==2], 
        {Tube[bonds[[#]], bondstyle[[4]], VertexColors->Table[{Specularity[White, 100], atomcol[[confP[[tuples[[#, ii]]]]]]}, {ii, 1, 2}]]&/@Range[Length[bonds]]}, 
        (*bonds in a single colour: *)
        True, {bondstyle[[3]], Tube[#, bondstyle[[4]]]&/@bonds}
    ]*)
    
    }, Boxed->False, SphericalRegion->True, ImageSize->500, Lighting->"Accent", 
    PlotLabel->Style[label, 24, Bold, Black]
    ]
    ];


    (*LatticeVectors, SymEquiv, SymAllAtoms, ImportCIF 
    taken from https://library.wolfram.com/infocenter/MathSource/9373/*)
    (*construct lattice vectors from lengths and angles*)
    LatticeVectors[lengths_,angles_]:=
    Module[{a,b,c,\[Alpha],\[Beta],\[Gamma],veca,vecb,vecc,c2,c3},
        {a,b,c}=N[lengths];
        {\[Alpha],\[Beta],\[Gamma]}=N[angles]/180*\[Pi];
        veca = a*{1,0,0};
        vecb = b*{Cos[\[Gamma]],Sin[\[Gamma]],0};
        vecc = c*{Cos[\[Beta]],c2,c3};
        c2 = (Cos[\[Alpha]]-Cos[\[Gamma]]*Cos[\[Beta]])/Sin[\[Gamma]];
        c3 = Sqrt[1-(Cos[\[Beta]])^2-(c2)^2];
        Chop[{veca,vecb,vecc}]
    ];


    (*find all atoms symmetrically equivalent to one specific atom*)
    SymEquiv[sym_,coord_,conf_,index_]:=
        Module[{xtemp,ytemp,ztemp,symmetries,unique,symout},
            (*apply all symmetries found for the system to the atom, returning all symmetrically equivalent atoms*)
            symmetries=sym/.{ImportCIF`Private`x->coord[[index,1]],ImportCIF`Private`y->coord[[index,2]],ImportCIF`Private`z->coord[[index,3]]};
            (*delete duplicate atoms*)
            unique = Round[#,.000001]-Floor[Round[#,.000001]]&/@symmetries;
            unique = DeleteDuplicates[unique,(Abs[#1[[1]]-#2[[1]]]<.0001&&Abs[#1[[2]]-#2[[2]]]<.0001&&Abs[#1[[3]]-#2[[3]]]<.0001)&];
            (*return atoms and their types*)
            symout = {#,conf[[index]]}&/@unique;
            symout
        ];


    (*find all unique atoms from the pool of symmetrical equivalents*)
    SymAllAtoms[sym_,coord_,conf_]:=
        Module[{symequiv},
            (*run SymEquiv over all values of 'index', then delete duplicate atoms*)
            symequiv = Flatten[SymEquiv[sym,coord,conf,#]&/@Range[Length[conf]],1];
            DeleteDuplicates[symequiv,(Abs[#1[[1,1]]-#2[[1,1]]]<.0001&&Abs[#1[[1,2]]-#2[[1,2]]]<.0001&&Abs[#1[[1,3]]-#2[[1,3]]]<.0001&&TrueQ[#1[[2]]==#2[[2]]])&]
    ];


    (*read a CIF file: *)
    ImportCIF[file_String/;(FileExistsQ[file] && TrueQ[FileExtension[file]=="cif"])]:=
    Module[{crystal,data,cifMol,unstring,lengths,angles,
            lattvec,coord,coordall,confall,
            confnonchem,chemnonchem,conf,knowntypes,chem,
            symopstrings,splitstrings,sym,symequiv,
            explicitInput,fullinput,symequivnonchem,
            explicitInputNonChemfullinputnonchem,chosenProperties,
            coordnonchem,whichlabel,label,atomKeys},
    chosenProperties = {EntityProperty["Element", "AtomicRadius"], 
        EntityProperty["Element", "CovalentRadius"], 
        EntityProperty["Element", "ElectronConfiguration"], 
        EntityProperty["Element", "ElectronConfigurationString"], 
        EntityProperty["Element", "IconColor"], 
        EntityProperty["Element", "MostCommonOxidationStates"], 
        EntityProperty["Element", "Name"], 
        EntityProperty["Element", "RefractiveIndex"], 
        EntityProperty["Element", "SpaceGroup"], 
        EntityProperty["Element", "StableIsotopes"], 
        EntityProperty["Element", "VanDerWaalsRadius"]};

    crystal = <||>;
    (*import it as a mathematica-molecule*)
    cifMol = Import[file, "Molecule"];
    AssociateTo[crystal, "molecule" -> cifMol];
    cifMol = List @@ cifMol;
    AssociateTo[crystal,"atoms" -> ((# -> <|"entity" -> Entity["Element", #]|>) & /@ (#[[1]] & /@ cifMol[[1]]))];
    crystal["atoms"] = Association[crystal["atoms"]];
    atomKeys = Keys[crystal["atoms"]];

    Do[(key = atomKeys[[i]];
        crystal["atoms"][key] = (Merge[{Association[#[[1]] -> #[[2]] & /@ 
                 Transpose[{(#[[2]] & /@ chosenProperties), 
                     EntityValue[crystal["atoms"][key]["entity"], 
                        chosenProperties]}]], crystal["atoms"][key]}, Total])), {i, 1, Length[atomKeys]}];
    
    (*import it as an Association*)
    cifText = Import[file, "Text"] <> "\n#";
    data = Association[ImportString[cifText, "CIF"]];
    crystal["space_group_HM"] = data["_symmetry_space_group_name_H-M"];
    
    (*function to convert to an expression, and remove uncertainty digits: *)
    unstring[x_]:= If[MatchQ[x,_String],
                  ToExpression[StringSplit[x,{"(",")","[","]","{","}",","}][[1]]],x];

    (*lattice vectors*)
    lengths = {data["_cell_length_a"], data["_cell_length_b"], data["_cell_length_c"]};
    angles = {data["_cell_angle_alpha"], data["_cell_angle_beta"], data["_cell_angle_gamma"]};
    lattvec = LatticeVectors[unstring/@lengths, unstring/@angles];

    (*atoms and their types*)
    coord = Transpose[{unstring/@data["_atom_site_fract_x"],
                     unstring/@data["_atom_site_fract_y"],
                     unstring/@data["_atom_site_fract_z"]}];
    conf = data["_atom_site_label"];
    confnonchem = conf;
    chemnonchem = DeleteDuplicates[confnonchem];
    confnonchem = ReplacePart[confnonchem,Flatten@Table[#->ii&/@(Flatten@Position[confnonchem,DeleteDuplicates[confnonchem][[ii]]]),{ii,1,Length@DeleteDuplicates[confnonchem]}]];
    
    (*separate proper elements from appendices: *)
    knowntypes = {"Wat","Uuu","Uut","Uus","Uuq","Uup","Uuo","Uuh","Uub",
        "O-H","Zr","Zn","Yb","Xe","Tm","Tl","Ti","Th","Te","Tc","Tb",
        "Ta","Sr","Sn","Sm","Si","Sg","Se","Sc","Sb","Ru","Rn",
        "Rh","Rg","Rf","Re","Rb","Ra","Pu","Pt","Pr","Po","Pm",
        "Pd","Pb","Pa","Os","Np","No","Ni","Ne","Nd","Nb","Na",
        "Mt","Mo","Mn","Mg","Md","Lu","Lr","Li","La","Kr","Ir","In",
        "Hs","Ho","Hg","Hf","He","Ge","Gd","Ga","Fr","Fm","Fe",
        "Eu","Es","Er","Dy","Ds","Db","Cu","Cs","Cr","Co",
        "Cn","Cm","Cl","Cf","Ce","Cd","Ca","Br","Bk","Bi","Bh",
        "Be","Ba","Au","At","As","Ar","Am","Al","Ag","Ac","Y","W",
        "V","U","T","S","P","O","N","K","I","H","F","D","C","B",""};
    Table[(conf[[ii]] = Select[knowntypes, StringMatchQ[conf[[ii]], #<>"*"]&][[1]])
         ,{ii,Range[Length[conf]]}];
    chem = DeleteDuplicates[conf];
    conf = ReplacePart[conf, Flatten@Table[#->ii&/@(Flatten@Position[conf, DeleteDuplicates[conf][[ii]]]),
                    {ii,1,Length@DeleteDuplicates[conf]}]];

    (*substance name: *)
    whichlabel = Select[{"_chemical_name_mineral", "_amcsd_formula_title", "_pd_phase_name", "_chemical_formula_structural","_chemical_formula_iupac","_chemical_formula_sum"}, KeyExistsQ[data,#]&];
    label = If[TrueQ[whichlabel=={}],FileBaseName[file], data[whichlabel[[1]]]];

    (*asymmetric units: *)
    explicitInput=<|"lattice"->lattvec,
                    "atomcoords"->coord,
                    "atomtypes"->conf,
                    "chemical"->chem,
                    "name"->label,
                    "file"->AbsoluteFileName[file],
                    "comment"->"atoms of the asymmetric unit with chemical types"|>;

    explicitInputnonchem = <|"lattice"->lattvec,
                             "atomcoords"->coord,
                             "atomtypes"->confnonchem,
                             "chemical"->chemnonchem,
                             "name"->label,
                             "file"->AbsoluteFileName[file],
                             "comment"->"atoms of the asymmetric unit with NONchemical types"|>;

    (*read symmetries: *)
    symopstrings = data[Select[{"_space_group_symop_operation_xyz","_symmetry_equiv_pos_as_xyz"}
                     ,KeyExistsQ[data,#]&][[1]]];
    splitstrings = StringSplit[#,","]&/@symopstrings;
    splitstrings = StringReplace[#,{"x"->"ImportCIF`Private`x",
                                  "y"->"ImportCIF`Private`y",
                                  "z"->"ImportCIF`Private`z"}]& /@ splitstrings;
    sym = ToExpression[splitstrings];

    (*add symmetrically equivalent atoms for chemical structures: *)
    symequiv = SymAllAtoms[sym,coord,conf];
    coordall = symequiv[[All,1]];
    confall = symequiv[[All,2]];
    fullinput = <|"lattice"->lattvec,
                "atomcoords"->coordall,
                "atomtypes"->confall,
                "chemical"->chem,
                "name"->label,
                "file"->AbsoluteFileName[file],
                "comment"->"complete cell with chemical types"|>;
    
    (*add symmetrically equivalent atoms for non-chemical structures: *)
    If[TrueQ[chem==chemnonchem],
        fullinputnonchem=fullinput;,
        symequivnonchem=SymAllAtoms[sym,coord,confnonchem];
        coordnonchem=symequivnonchem[[All,1]];
        confnonchem=symequivnonchem[[All,2]];
        fullinputnonchem=<|"lattice"->lattvec,
                            "atomcoords"->coordnonchem,
                            "atomtypes"->confnonchem,
                            "chemical"->chemnonchem,
                            "name"->label,
                            "comment"->"complete cell with NONchemical types"|>;
        ];
    
    (*crystal["explitic_non_chemical"] = explicitInputnonchem;
    crystal["full_input_non_chemical"] = fullinputnonchem;
    crystal["explicit_input"] = explicitInput;*)
    crystal["crystal_structure"] = fullinput;
    crystal["source_file"] = FileNameSplit[AbsoluteFileName[file]][[-1]];
    cystal["data"] = data;
    Dataset[crystal]
    ];
    
AutoRotate[gr_Graphics3D,rate_:5]:=DynamicModule[{vp,va,vv,vc},
    {vp,va,vv,vc}=gr~AbsoluteOptions~#~OptionValue~#&@{ViewPoint,ViewAngle,ViewVertical,ViewCenter};
    (*from https://mathematica.stackexchange.com/questions/39398/how-to-make-a-3d-plot-auto-rotate*)
    Overlay[{
        Show[gr,
            SphericalRegion->True,
            ViewPoint->Dynamic[RotationMatrix[Clock[2 \[Pi],rate],vv].vp],
            ViewAngle->Dynamic[va],
            ViewVertical->Dynamic[vv],
            ViewCenter->Dynamic[vc]
        ]},
    All,1]];
    
End[] (*End private*)
EndPackage[]
import Foundation
import MathTools

// this handles special case of a membrane

public class UnitCell {

    var origin:Vector
    var dimensions:[Vector]
    var deltas:[Double]
    var buffer:[Double]
    var units:[Vector]
    var sizes:[Double]
    var membraneaxis:AXES
    var imagecoords:[Vector]?
    var imageindices:[Int:[Int]]?
    var inverseindices:[Int:Int]?


    public init(_ origin:Vector, _ dimensions:[Vector], _ buffer:[Double], _ deltas:[Double],  _ membraneaxis:AXES )  {
        self.origin = origin
        self.dimensions = dimensions
        self.deltas = deltas
        // note that unit vector is an optional
        self.units = dimensions .map { $0.unit()! }
        self.sizes = dimensions .map { $0.length() }
        // translate buffer distances into cell coordinates
        self.buffer = zip( buffer, self.sizes ) .map { $0 / $1 }
        self.membraneaxis = membraneaxis

        self.imagecoords = nil
        self.imageindices = nil
        self.inverseindices = nil

    }

    public func cellcoords( _ coord:Vector ) -> [Double] {
        
        let disp = coord.sub(self.origin)
        let cellc = (0..<3) .map { disp.dot(self.units[$0]) / self.sizes[$0] }

        return cellc
    }

    public func packcell( _ atompos:[Vector] ) -> [Vector] {

    var transPos = [Vector]()

    // this assumes at most one image of each atom is present

    for coord in atompos {

            let cellcoords = self.cellcoords(coord)
            var disp = coord.sub(self.origin)

            for ax in 0..<3 {
                
                if cellcoords[ax] < 0.0 {
                    disp = disp.add(self.dimensions[ax])
                }

                if cellcoords[ax] > 1.0 {
                    disp = disp.sub(self.dimensions[ax])
                }
            }

            transPos.append(self.origin.add(disp))

        }

        return transPos

    }

    // this adds reflections of atoms that are within 'buffer' cell coordinate of the cell boundaries
    // in the 'non-cell membrane' directions
    // 
    // a coordinate can have one image if it is within one buffer boundary, or three if it is close to two boundaries

    public func buildimages( _ coords:[Vector] ) -> ([Vector], [Int:[Int]], [Int:Int]) {

        var bcoords = [Vector]()
        var bmap = [Int:[Int]]()
        var binvmap = [Int:Int]()

        

        for (aidx,coord) in coords.enumerated() {
            let cellc = self.cellcoords( coord )
            var disps = [Vector]()
    

            for ax in 0..<3 {
                if ax == self.membraneaxis.rawValue {
                    continue
                }

                if cellc[ax] < self.buffer[ax] {
                    disps.append(self.dimensions[ax].scale(1.0)) 

                }

                if 1.0 - cellc[ax] < self.buffer[ax] {
                    disps.append(self.dimensions[ax].scale(-1.0))
                }
            }

            if disps.count > 0  {
                if bmap[aidx] == nil {
                    bmap[aidx] = [Int]()
                }
                // apply displacements single and in combination
                var apply = [Vector]()
                if disps.count == 2 {
                    apply.append(disps[0])
                    apply.append(disps[1])
                    apply.append(disps[0].add(disps[1]))
                }
                else {
                    apply.append(disps[0])
                }

                for disp in apply {
                    bmap[aidx]!.append(bcoords.count)
                    bcoords.append(coord.add(disp))
                }
                
            }
        }

        // for convenience make inverse map (image index to original atom index)

        for aidx in bmap.keys {
            for imgidx in bmap[aidx]! {
                binvmap[imgidx] = aidx
            }
        }
        return (bcoords, bmap, binvmap)

    }

    public func setimages( _ imagecoords:[Vector], _ imageindices:[Int:[Int]], _ inverseindices:[Int:Int]) {
        self.imagecoords = imagecoords
        self.imageindices = imageindices
        self.inverseindices = inverseindices
    }

}



public func membraneCoordinates(_ incoordinates:[Vector], _ inradii:[Double], _ proberad:Double, _ unitcell:UnitCell ) 
    -> ([Vector], ([Vector], [Int:[Int]], [Int:Int]), [Double])  {

        var packedcoords = unitcell.packcell( incoordinates )

        let imgdata = unitcell.buildimages(packedcoords) 

        unitcell.setimages( imgdata.0, imgdata.1, imgdata.2 )

        var imgwithidx = [(Int,Int)]()

        for (aidx,imgindices) in unitcell.imageindices! {
            for imgidx in imgindices {
                imgwithidx.append( (imgidx,aidx ))
            }
        }

        let imginorder = imgwithidx .sorted { $0.0 < $1.0 }

        let imgradii = imginorder .map { inradii[$0.1] }

        return (packedcoords, imgdata, imgradii)

}

public func processMembraneProbes( _ probes:[Probe], _ proberad:Double, _ unitcell:UnitCell ) 
            -> [Probe] {

    // ** in first step, remove probes that lie outside the unit cell ** 
    // I think this is counterproductive - I am potentially losing information
    // from probes that span one unit cell to the next 
    //
    // instead of building 'image probes', just retain all probes within buffer of both boundaries at the outset
    // I may have been worrying about 'double counting' density, but in latest rev that is not possible

    var keepprobes = [Probe]()


    for probe in probes {
        let cellc = unitcell.cellcoords(probe.center)
        var outside = false
        for ax in 0..<3 {
            if ax == unitcell.membraneaxis.rawValue {
                continue
            }
            if cellc[ax] > 1.0 + unitcell.buffer[ax] || cellc[ax] < 0.0 - unitcell.buffer[ax] {
                outside = true
            }
        }
        if outside {
            continue
        }

        keepprobes.append(probe)

    }

    // second step is to generate image probes using the unit cell buffer dimensions
    // do not worry about duplicates at this point

    /*
    let probecoords = keepprobes .map { $0.center }

    let probeimgdata = unitcell.buildimages( probecoords )

    let probeinvimg = probeimgdata.2 

    var bufferPROBES = [Probe]() 

    for imgidx in probeinvimg.keys {
        let probidx = probeinvimg[imgidx]!
        let probe = keepprobes[probidx]
        let bufprobe = Probe(center:probeimgdata.0[imgidx], radius:probe.proberadius, atoms:probe.atoms, 
            singleton:probe.singleton, clockwise:probe.clockwisecontour)
        bufferPROBES.append(bufprobe)


    }
    
    let totPROBES = keepprobes + bufferPROBES
    */

    return keepprobes


}



public func processMembraneTri( _ VERTICES:[Vector], _ NORMALS:[Vector], _ FACES:[[Int]], _ unitcell:UnitCell ) 
        -> ([Vector], [Vector], [[Int]]) {

    var inelems = [Int]()
    var outelems = [Int]()
    var crosselems = [Int]()

    let vertexcellcoords = VERTICES .map { unitcell.cellcoords($0) }

    for (fidx,face) in FACES.enumerated() {
        var inside = 0
        var outside = 0
        for vertex in face .map ({ vertexcellcoords[$0] }) {
            var out = false

            for ax in 0..<3 {
                if ax == unitcell.membraneaxis.rawValue { continue }
                if vertex[ax] < 0.0 || vertex[ax] > 1.0 {
                    out = true
                    break
                }
            }
            
            if out {
                outside += 1
            }
            else {
                inside += 1
            }

        } 
        if outside == 0 {
            inelems.append(fidx)
        }
        else if inside == 0 {
            outelems.append(fidx)
        }
        else if outside > 0 && inside > 0 {
            crosselems.append(fidx)
        }
    } 

    // not sure that initial approach of simply deleting elements at one side is adequate - what happens upon decimation??
    // 
    // instead add edges to make boundary elements conform to unit cell boundaries; means adding vertices, edges and faces

    // for oriented element [0, 1, 2, 0] :
    // if two vertices separated by u.c. boundary, insert a vertex between 
    // e.g. if 2 and 0 on opp sides, have [0, 1, 2, 20, 0 ]
    // have faces [0, 1, 2] + [0, 2, 20]
    //
    // if 0, 2 opp 1 : [0, 01, 1, 2, 20, 0], elements [0, 01, 20], [01, 1, 20], [1, 2, 20]
    //
    // special condition - does face contain unit cell 'corner' ? This gets really complicated; maybe see first if 
    // smoothing and decimation looks OK with the ragged edge.

    /*
    var newfaces = [[Int]]()

    for fidx in crosselems {

        let insideverts =  FACES[fidx] .filter { vertexcellcoords[$0] <= 1.0 || vertexcellcoords[$0] >= 0.0 }
        let outsideverts = FACES[fidx] .filter { vertexcellcoords[$0] > 1.0 || vertexcellcoords[$0] < 0.0 }

        var extendedface = [Int]()

        for i in 0..<3 {
            let j = i < 2 ? i + 1 : 0

            let vi = FACES[fidx][i]
            let vj = FACES[fidx][j]

            if (insideverts.contains(vi) && outsideverts.contains(vj)) || (insideverts.contains(vj) && outsideverts.contains(vi)) {
                
            }
        }



    } 

    */



    // delete crossover elements at 'right' or 'top' of unit cell --- 
    // in brief, keep any element with all vertex coordinates less than 1.0

    var keepelems = Array(inelems) 


    for fidx in crosselems {
        // just keep for now, sort out after I see the results
        keepelems.append(fidx)
        /*
        var keep = true 
        let vertcoord = FACES[fidx] .map { vertexcellcoords[$0] }
        
        for (aidx,coord) in vertcoord.enumerated() {
            if aidx == unitcell.membraneaxis.rawValue { continue }
            if coord[aidx] > 1.0 {
                keep = false
                break
            }
        }

        if keep {
            keepelems.append(fidx)
        }
        */
    }

    var vertexindices = [Int]()

    for fidx in keepelems {
        vertexindices += FACES[fidx]
    }

    let keepvertexset = Set(vertexindices)

    var keepvertexindices = Array(keepvertexset) 
    keepvertexindices = keepvertexindices .sorted { $0 < $1 }

    // translate old to new vertex indices 

    var translate = Array( repeating:-1 , count:VERTICES.count )

    for (newidx,oldidx) in keepvertexindices.enumerated() {
        translate[oldidx] = newidx
    }

    // vertex data for export 

    let exportVertices = keepvertexindices .map { VERTICES[$0] }
    let exportNormals =  keepvertexindices .map { NORMALS[$0] }

    var exportElems = [[Int]]()

    for fidx in keepelems {
        let v0 = FACES[fidx][0]
        let v1 = FACES[fidx][1]
        let v2 = FACES[fidx][2]
        exportElems.append( [translate[v0], translate[v1], translate[v2]])
    }

    return (exportVertices, exportNormals, exportElems)


}



public func mergeSurfaceComponents( _ vertA:[Vector], _ vertB:[Vector], _ normA:[Vector], _ normB:[Vector],
    _ faceA:[[Int]], _ faceB:[[Int]]) -> ([Vector],[Vector],[[Int]]) {

        var outverts = [Vector]() 

        var outnorms = [Vector]()

        outverts += vertA 
        outverts += vertB

        outnorms += normA 
        outnorms += normB 

        let vertexmapB = (0..<vertB.count) .map { $0 + vertA.count }

        var outfaces = [[Int]]()

        outfaces += faceA 

        outfaces += faceB .map { [vertexmapB[$0[0]], vertexmapB[$0[1]], vertexmapB[$0[2]]] }

        return ( outverts, outnorms, outfaces ) 
    } 



public func membraneSurfaceComponents( _ SUBVERTICES:[[Vector]], _ SUBNORMALS:[[Vector]], _ SUBFACES:[[[Int]]], 
            _ unitcell:UnitCell)
             ->  ([[Vector]], [[Vector]], [[[Int]]])  {

        // strategy - take random sample of biggest components, for a sample of elements get normal along membrane axis, 
        // average coordinate in direction of axis; expect highest and lowest with opposite normals, inner two with likwise 
        // opposite normals ; merge outer and inner into single components ('0' is outer, '1' is inner)

        var OUTVERTICES = [[Vector]]()
        var OUTNORMALS = [[Vector]]()
        var OUTFACES = [[[Int]]]()

        var maxPositions = [Double]()

        let ax = unitcell.membraneaxis.rawValue
        

        for isurf in 0..<4 {
            

            var axcoors = SUBVERTICES[isurf] .map { $0.coords[ax] }
            
            maxPositions.append( axcoors.max()! )

        }

        

        var compAxisPos = (0..<4) .map { (maxPositions[$0],$0) } 

        compAxisPos = compAxisPos .sorted { $0.0 < $1.0 }

        let bottom0 = compAxisPos[0].1
        let bottom1 = compAxisPos[1].1
        let top1 = compAxisPos[2].1
        let top0 = compAxisPos[3].1

        
        // membrane components are sorted - merge top.0 and bottom.0, top.1 and bottom.1 ; any
        // other components are copied

        let membrane0 = mergeSurfaceComponents( SUBVERTICES[top0], SUBVERTICES[bottom0], 
                SUBNORMALS[top0], SUBNORMALS[bottom0], SUBFACES[top0], SUBFACES[bottom0] )

        let membrane1 = mergeSurfaceComponents( SUBVERTICES[top1], SUBVERTICES[bottom1], 
                SUBNORMALS[top1], SUBNORMALS[bottom1], SUBFACES[top1], SUBFACES[bottom1] )

        OUTVERTICES.append(membrane0.0)
        OUTVERTICES.append(membrane1.0)

        OUTNORMALS.append(membrane0.1)
        OUTNORMALS.append(membrane1.1)

        OUTFACES.append(membrane0.2)
        OUTFACES.append(membrane1.2)

        if SUBVERTICES.count > 4 {
            for isurf in 4..<SUBVERTICES.count {
                OUTVERTICES.append(SUBVERTICES[isurf])
                OUTNORMALS.append(SUBNORMALS[isurf])
                OUTFACES.append(SUBFACES[isurf])
            }
        }

        return (OUTVERTICES,OUTNORMALS,OUTFACES)


    }
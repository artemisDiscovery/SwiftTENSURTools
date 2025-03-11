import Foundation
import MathTools

// this handles special case of a membrane

public class UnitCell {

    var origin:Vector
    var dimensions:[Vector]
    var levelspacings:[Double]
    var griddeltas:[Double]
    var buffer:[Double]
    var units:[Vector]
    var sizes:[Double]
    var membraneaxis:AXES
    var imagecoords:[Vector]?
    var imageindices:[Int:[Int]]?
    var inverseindices:[Int:Int]?


    public init( origin:Vector, dimensions:[Vector], buffer:[Double], levelspacings:[Double], griddeltas:[Double],   membraneaxis:AXES )  {
        self.origin = origin
        self.dimensions = dimensions
        self.griddeltas = griddeltas
        self.levelspacings = levelspacings
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

    public func inside( _ coords:[Vector], tol:Double=0.0001 ) -> Bool {
        let cellc = self.cellcoords( coords )

        for ax in 0..<3 {
            if cellc[ax] < -tol || cellc[ax] > 1.0 + tol {
                return false
            }
        }

        return true
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
            -> ([Probe],[Probe],[Probe]) {

    // ** in first step, remove probes that lie outside the unit cell ** 
    // I think this is counterproductive - I am potentially losing information
    // from probes that span one unit cell to the next 
    //
    // instead of building 'image probes', just retain all probes within buffer of both boundaries at the outset
    // I may have been worrying about 'double counting' density, but in latest rev that is not possible

    var keepprobes = [Probe]()
    var bufferPROBES = [Probe]() 
    var totPROBES = [Probe]() 


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

    return (keepprobes, bufferPROBES, totPROBES)


}



public func processMembraneTri( VERTICES:[Vector], NORMALS:[Vector], FACES:[[Int]], PROBES:[Probe], unitcell:UnitCell  ) 
        -> (vertices:[Vector], normals:[Vector], faces:[[Int]], surfacetype:SurfaceType) {


    // BEFORE doing the cut, determine connected components; 

    let time0 = Date().timeIntervalSince1970
    let components = surfaceComponentsDFS( faces:FACES, numvertices:VERTICES.count )
    let time1 = Date().timeIntervalSince1970

    print("\ntime to find connected comoponents by DFS = \(time1 - time0)")

    // assume that four biggest (first) components are the membrane probe-centered and reentrant
    //
    // we want at this point only the merged reentrant, and with all elements retained of the corresponding
    // components

    var topPCIdx:Int?
    var bottomPCIdx:Int?
    var topReentIdx:Int?
    var bottomReentIdx:Int?

    var minMaxCompZ = [(min:Double,max:Double)]()

    for component in largestOpen.enumerated() {
        let compZ = component.vertices .map { $0.coords[2] }
        let minZ = compZ.min()!
        let maxZ = compZ.max()!
        minMaxCompZ.append((min:minZ, max:maxZ))
    }

    for (cidx,zdata) in minMaxCompZ.enumerated() {
        if zdata.max > probeZMax {
            topPCIdx = cidx
        }
        else if zdata.min < probeZMin {
            bottomPCIdx = cidx
        }
    }

    if topPCIdx == nil || bottomPCIdx == nil {
        print("\nmembrane surface components FAILURE, cannot identify both top and bottom probe-centered membrane components" )
        return nil 
    }

    var reentrantComponents = [] 

    for cidx in 0..<4 {
        if ![topPCIdx!, bottomPCIdx!].contains(cidx) {
            reentrantComponents.append(cidx)
        }
    }

    for order in [(0,1),(1,0)] {
        let top = reentrantComponents[order[0]]
        let bottom = reentrantComponents[order[1]]
        if minMaxCompZ[top].min > minMaxCompZ[bottom].max {
            topReentIdx = top
            bottomReentIdx = bottom
            break
        }
    }

    if topReentIdx == nil || bottomReentIdx == nil {
        print("\nmembrane surface components FAILURE, cannot identify both top and bottom reentrant membrane components" )
        return nil 
    }

    var vertexToFaces = Array( repeating:[Int](), count:VERTICES.count )

    for (fidx,face) in FACES.enumerated() {
        face .map { vertexToFaces[$0].append(fidx) }
    }

    var membranevertices = Array( repeating:false, count:VERTICES.count )

    for cidx in reentrantComponents {

        _ = components[cidx] .map { membranevertices[$0] = true }
    }


    // only keep membrane faces/vertices inside the unit cell 

    let insidevertices = VERTICES .map { unitcell.inside( $0 ) }

    let keepvertices = zip(insidevertices, membranevertices) .map { $0 && $1 }

    var membranefaces = Array( repeating:true, count:FACES.count )

    _ = FACES.enumerated() .map { ($0) in 
        if !keepvertices[$0.element[0]] || !keepvertices[$0.element[1]] || !keepvertices[$0.element[2]] {
            membranefaces[$0.offset] = false 
        }
    }

    // renumber to return 

    let exportvertexindices = 0..<VERTICES.count .filter { keepvertices[$0] }

    var oldToNewIndex = Array(repeating:-1, count:VERTICES.count )

    for (newidx,oldidx) in exportvertexindices.enumerated() {
        oldToNewIndex[oldidx] = newidx 
    } 

    let exportvertices = VERTICES.enumerated() .filter { keepvertices[$0.offset] } .map { $0.element } 
    let exportnormals  =  NORMALS.enumerated() .filter { keepvertices[$0.offset] } .map { $0.element }

    let exportfaces = FACES.enumerated() .filter { membranefaces[$0.offset] } .map { $0.element .map { oldToNewIndex[$0]} }


 

    return (vertices:exportvertices, normals:exportnormals, faces:exportfaces, surfacetype:SurfaceType.undeterminedOpen)


}






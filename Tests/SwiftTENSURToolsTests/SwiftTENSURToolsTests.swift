
import XCTest
@testable import SwiftTENSURTools

import Foundation 
import MathTools

enum tensurError: Error {
    case testImportError
    case parseError
}

func readCoordinateFile( _ path:String ) throws -> [Vector] {

    var lines:[String]?

    do {
        lines = try String(contentsOfFile:path).split { $0 .isNewline } .map { String($0) }
    }
    catch {
        print("exception reading data file")
        throw tensurError.testImportError
    }
    
    print("in readCoordinateFile, have \(lines!.count) lines imported from \(path)")

    var coords = [Vector]()

    for lidx in 0..<lines!.count {
        let tokens = lines![lidx].split(separator: " ")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
        
        if tokens.count < 3 {
            continue
        }

        let x = try Double(tokens[0])
        let y = try Double(tokens[1])
        let z = try Double(tokens[2])

        if x == nil || y == nil || z == nil {
                throw tensurError.parseError
        }

        let coord = Vector([x!,y!,z!])

        coords.append(coord)

    }

    return coords
}

func importOBJFile( _ path:String ) throws -> ([Vector],[Vector],[[Int]]) {

    var lines:[String]?

    do {
        lines = try String(contentsOfFile:path).split { $0 .isNewline } .map { String($0) }
    }
    catch {
        print("exception reading data file")
        throw tensurError.testImportError
    }

    var VERTICES = [Vector]()
    var NORMALS = [Vector]()
    var FACES = [[Int]]()

    for lidx in 0..<lines!.count {

        let tokens = lines![lidx].split(separator: " ")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
        
        if tokens.count < 1 {
            continue
        }

        var linetype:String?

        linetype = try String(tokens[0])

        if linetype == nil {
            throw tensurError.parseError
        }

        if linetype == "v" || linetype == "vn" {
            var coords = [Double]()
            var nextcoord:Double?

            for idx in 0..<3 {
                nextcoord = try Double(tokens[1 + idx])
                if nextcoord == nil {
                    throw tensurError.parseError
                }
                coords.append(nextcoord!)
            }
            if linetype == "v" {
                VERTICES.append(Vector(coords))
            }
            else {
                NORMALS.append(Vector(coords))
            }
            
        }
        else if linetype == "f" {
            var indices = [Int]()
            var nextidx:Int?

            for idx in 0..<3 {
                nextidx = try Int(tokens[1 + idx])
                if nextidx == nil {
                    throw tensurError.parseError
                }
                indices.append(nextidx! - 1)
            }

            FACES.append(indices)


        }

    }

    return (VERTICES,NORMALS,FACES)

}
func readRadiiFile( _ path:String ) throws -> [Double] {

    var lines:[String]?

    do {
        lines = try String(contentsOfFile:path).split { $0 .isNewline } .map { String($0) }
    }
    catch {
        print("exception reading data file")
        throw tensurError.testImportError
    }
    

    var radii = [Double]()

    for lidx in 0..<lines!.count {
        let tokens = lines![lidx].split(separator: "\t")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
        
        if tokens.count < 1 {
            continue
        }

        let rad = try Double(tokens[0])

        if rad == nil  {
                throw tensurError.parseError
        }

        radii.append(rad!)

    }

    return radii
}


func readCoordinates( _ lines:[String], _ index : inout Int ) throws -> [Vector] {
    let ncoord = Int(lines[index])!

    var coords = [Vector]()


        for _ in 0..<ncoord {
            index += 1
            let tokens = lines[index].split(separator: "\t") 
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            //print("tokens : \(tokens)")
            
            let x = try Double(tokens[0])
            let y = try Double(tokens[1])
            let z = try Double(tokens[2])

            if x == nil || y == nil || z == nil {
                throw tensurError.parseError
            }

            let coord = Vector([x!,y!,z!])
            coords.append(coord)
            
            
        }
    
    // get to next line 
    index += 1
    return coords 
}

func readRadii( _ lines:[String], _ index : inout Int ) -> [Double] {
    let nrad = Int(lines[index])!

    var radii = [Double]()

    for _ in 0..<nrad {
        index += 1
        let rad = Double(lines[index].trimmingCharacters(in: .whitespacesAndNewlines))!
        radii.append(rad)
    }
    // get to next line 
    index += 1
    return radii 
}

func readCircles(_ lines:[String], _ index : inout Int) -> ([Double],[[AtomCircle]]) {
    let tokens = lines[index].split(separator:"\t")
    let nlevels = Int(tokens[0])!

    var levels = [Double]()
    var circlesForLevels = [[AtomCircle]]()

    for _ in 0..<nlevels {
        index += 1
        let tokens = lines[index].split(separator:"\t")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
        let level = Double(tokens[0])!
        levels.append(level)
        var circles = [AtomCircle]()
        let ncircles = Int(tokens[1])!

        for _ in 0..<ncircles {
            index += 1
            var tokens = lines[index].split(separator:"\t")
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            let atom = Int(tokens[0])!
            index += 1
            tokens = lines[index].split(separator:"\t")
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            let x = Double(tokens[0])!
            let y = Double(tokens[1])!
            let z = Double(tokens[2])!
            let r = Double(tokens[3])!
            
            circles.append(AtomCircle(atom, Vector([x,y,z]), r, AXES.X))
        }

        circlesForLevels.append(circles)
    }

    index += 1 

    return (levels,circlesForLevels)

}

func readContours(_ lines:[String], _ index : inout Int) -> [[[(Int,Vector?,Int?,Vector?,Int?)]]] {
    let tokens = lines[index].split(separator:"\t")
        .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
    let nlevels = Int(tokens[0])!

    var contoursForLevels = [[[(Int,Vector?,Int?,Vector?,Int?)]]]()

    for _ in 0..<nlevels {
        index += 1
        var tokens = lines[index].split(separator:"\t")
            .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
        
        // contour = [(circleidx,Vector?,Int?,Vector?,Int?)]
        var contours = [[(Int,Vector?,Int?,Vector?,Int?)]]()
        let ncontours = Int(tokens[1])!

        for _ in 0..<ncontours {
            index += 1
            tokens = lines[index].split(separator:"\t")
                .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
            let narcs = Int(tokens[0])!
            var contour = [(Int,Vector?,Int?,Vector?,Int?)]()
            for _ in 0..<narcs {
                index += 1
                tokens = lines[index].split(separator:"\t")
                    .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
                let cidx = Int(tokens[0])!
                if lines[index].contains("nonsingleton") {
                    index += 1
                    tokens = lines[index].split(separator:"\t")
                        .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
                    let xI = Double(tokens[0])!
                    let yI = Double(tokens[1])!
                    let zI = Double(tokens[2])!
                    let circleI = Int(tokens[3])!

                    index += 1
                    tokens = lines[index].split(separator:"\t")
                        .map { $0.trimmingCharacters(in: .whitespacesAndNewlines)}
                    let xJ = Double(tokens[0])!
                    let yJ = Double(tokens[1])!
                    let zJ = Double(tokens[2])!
                    let circleJ = Int(tokens[3])!

                    contour.append( (cidx,Vector([xI,yI,zI]),circleI,Vector([xJ,yJ,zJ]),circleJ) )
                }
                else {
                    contour.append( (cidx,nil,nil,nil,nil) )
                }
            }
            contours.append(contour)
        }

        contoursForLevels.append(contours)
    }

    index += 1 

    return contoursForLevels 

}


func importData(_ path:String) throws -> 
    ([Vector]?,[Double]?,Double?,[Double]?,[[AtomCircle]]?,[[[(Int,Vector?,Int?,Vector?,Int?)]]]?) {

    var input:[String]?

    do {
        input = try String(contentsOfFile:path).split { $0 .isNewline } .map { String($0) }
    }
    catch {
        print("exception reading data file")
        throw tensurError.testImportError
    }
    
    let lines = input!


    var lidx = 0 
    var coordinates:[Vector]?
    var radii:[Double]?
    var probeRad:Double?
    var circlesForLevels:[[AtomCircle]]?
    var levels:[Double]?

    var contoursForLevels:[[[(Int,Vector?,Int?,Vector?,Int?)]]]?

    while lidx < lines.count {
        if lines[lidx].contains("# coordinates") {
            lidx += 1
            do {
                coordinates = try readCoordinates(lines, &lidx )
            }
            catch {
                print("exception readCoordinates")
            }
            
            continue
        }

        if lines[lidx].contains("# radii") {
            lidx += 1
            radii = readRadii(lines, &lidx )
            continue
        }

        if lines[lidx].contains("# probeRad") {
            lidx += 1
            probeRad = Double(lines[lidx].trimmingCharacters(in: .whitespacesAndNewlines))
            lidx += 1
            continue
        }

        if lines[lidx].contains("# atom circles") {
            lidx += 1
            let cdata = readCircles(lines, &lidx )
            levels = cdata.0 
            circlesForLevels = cdata.1
            continue
        }

        if lines[lidx].contains("# contours") {
            lidx += 1
            contoursForLevels = readContours(lines, &lidx )
            
            continue
        }

        lidx += 1


    }

    return (coordinates,radii,probeRad,levels,circlesForLevels,contoursForLevels)


}

final class SwiftTENSURToolsTests: XCTestCase {
    func testExample() throws {
        // XCTest Documentation
        // https://developer.apple.com/documentation/xctest

        // Defining Test Cases and Test Methods
        // https://developer.apple.com/documentation/xctest/defining_test_cases_and_test_methods
    }

    // disable test by simply taking 'test' off the front of the function name !
    // https://pacugindre.medium.com/disable-xctests-on-swift-package-manager-quick-and-dirty-c69e94d8a2f0


    func testMembrane() throws {

        let buffer = [4.0,4.0,4.0]
        var unitcell = UnitCell(Vector([0.0,0.0,0.0]), 
            [Vector([61.864 , 0.0 , 0.0]), Vector([0.0 , 61.453 , 0.0]), Vector([0.0 , 0.0 , 96.120])], 
            buffer, AXES.Z ) 

        // import the mesh we saved in another test, here we just test the decompostion into components 

        var surfdata:([Vector],[Vector],[[Int]])?

        do {
            surfdata = try importOBJFile("membrane_proc.obj")
        }
        catch {
            print("exception importing object file ")
            throw tensurError.testImportError
        }

        let VERTICES = surfdata!.0
        let NORMALS = surfdata!.1
        let FACES = surfdata!.2

        print("\nimported surface has \(VERTICES.count) vertices , \(FACES.count) faces")
        
        // copy this from SwiftTENSUR, identifies subsurface components

        var adjacency = [Set<Int>]()

        for _ in 0..<VERTICES.count {
            adjacency.append(Set<Int>())
        }

        for f in FACES {
            adjacency[f[0]].insert(f[1])
            adjacency[f[1]].insert(f[0])
            adjacency[f[0]].insert(f[2])
            adjacency[f[2]].insert(f[0])
            adjacency[f[1]].insert(f[2])
            adjacency[f[2]].insert(f[1])

        }
        var STACK = [[Int]]()

        var visited = Array(repeating:false, count:VERTICES.count)
        var component = Array(repeating:-1, count:VERTICES.count)

        var currentComponent = -1
        var unassigned:Int?

        while true {
            //find first unassigned vertex

            unassigned = nil

            for iv in 0..<VERTICES.count {
                if component[iv] < 0 {
                    unassigned = iv
                    break
                }
            }
            if unassigned == nil {
                break
            }

            currentComponent += 1

            STACK.append([unassigned!,currentComponent])

            while STACK.count > 0 {
                let data = STACK.popLast()!
                if !visited[data[0]] {
                    visited[data[0]] = true 
                    component[data[0]] = data[1]
                    for iv in adjacency[data[0]] {
                        STACK.append([iv,currentComponent])
                    }
                }
            }
        }

        var components = [[Int]]()

        for c in 0..<(currentComponent+1) {
            components.append([Int]())
        }

        for (iv,c) in component.enumerated() {
            components[c].append(iv)
        }


        print("\nsurface has \(components.count) components")

        var SUBVERTICES = [[Vector]]()
        var SUBNORMALS = [[Vector]]()
        var SUBFACES = [[[Int]]]()

        // sort components by decreasing size

        components = components .sorted { $0.count > $1.count }

        for comp in components {
            let subvertindices = comp .sorted { $0 < $1 }
            var vertexmap = Array(repeating:-1, count:VERTICES.count)

            _ = subvertindices.enumerated() .map { vertexmap[$0.1] = $0.0 }

            let subvertices = subvertindices .map { VERTICES[$0] }
            let subnormals = subvertindices .map { NORMALS[$0] }

            let subfaces = FACES .filter { vertexmap[$0[0]] >= 0 } 
                .map { [vertexmap[$0[0]], vertexmap[$0[1]], vertexmap[$0[2]]] }

            SUBVERTICES.append( subvertices )
            SUBNORMALS.append( subnormals )
            SUBFACES.append( subfaces )
        }

        // merge membrane components

        let membranecomp = membraneSurfaceComponents( SUBVERTICES, SUBNORMALS, SUBFACES, unitcell )

        // component 0 is 'outside', component 1 is part we want 

        let VERTICES0 = membranecomp.0[0]
        let NORMALS0 = membranecomp.1[0]
        let FACES0 = membranecomp.2[0]
        let VERTICES1 = membranecomp.0[1]
        let NORMALS1 = membranecomp.1[1]
        let FACES1 = membranecomp.2[1]

        var url = URL(fileURLWithPath: "./membrane_0.obj")

        var outstr = ""

        for vertex in VERTICES0 {
            outstr += "v \(vertex.coords[0]) \(vertex.coords[1]) \(vertex.coords[2])\n"
        }

        for normal in NORMALS0 {
            outstr += "vn \(normal.coords[0]) \(normal.coords[1]) \(normal.coords[2])\n"
        }
        
        for face in FACES0 {
            outstr += "f \(face[0]+1) \(face[1]+1) \(face[2]+1)\n"
        }
        

        do {
            try outstr.write(to: url, atomically: true, encoding: String.Encoding.utf8)
        } catch {
            print("error writing file \(url)")
        }

        url = URL(fileURLWithPath: "./membrane_1.obj")

        outstr = ""

        for vertex in VERTICES1 {
            outstr += "v \(vertex.coords[0]) \(vertex.coords[1]) \(vertex.coords[2])\n"
        }

        for normal in NORMALS1 {
            outstr += "vn \(normal.coords[0]) \(normal.coords[1]) \(normal.coords[2])\n"
        }
        
        for face in FACES1 {
            outstr += "f \(face[0]+1) \(face[1]+1) \(face[2]+1)\n"
        }
        

        do {
            try outstr.write(to: url, atomically: true, encoding: String.Encoding.utf8)
        } catch {
            print("error writing file \(url)")
        }



}



 
    func donttestMembrane() throws {
        var coordinates:[Vector]?
        var radii:[Double]?

        do {
            coordinates = try readCoordinateFile( "Tests/DATA/coords.txt" )
        }
        catch {
            print("exception calling readCoordinateFile")
            throw tensurError.testImportError
        }

        do {
            radii = try readRadiiFile( "Tests/DATA/radii.txt" )
        }
        catch {
            print("exception calling readRadiiFile")
            throw tensurError.testImportError
        }

        // need unit cell : 61.864   61.453   96.120

        // unit cell is oriented in z-direction, make contours perpendicular to X and Y axes
        //
        // for X, the symmetry operation is 0, 61.864
        // for Y, 1, 61.453
        // make the buffer 4 (all radii < 2, proberad = 1.58)
        // I think the origin should be 0.,0.,0., as the coordinate means are about 1/2 the box dimensions
        // origin = 0., 0., 0.
        //
        // dimensions from max - min coordinate are quite different from the pdb unit cell, just reflects irregularity I think ? 
        // May need to make translations and check visually

        // 4841  C1  DPPC  404       6.559  31.299  64.890 
        // HETATM 4881  N41 DPPC  404       6.909  31.549  66.280 
        
        let buffer = [4.0,4.0,4.0]
        var unitcell = UnitCell(Vector([0.0,0.0,0.0]), 
            [Vector([61.864 , 0.0 , 0.0]), Vector([0.0 , 61.453 , 0.0]), Vector([0.0 , 0.0 , 96.120])], 
            buffer, AXES.Z ) 

        let probeRad = 1.58

        print("imported \(coordinates!.count) coordinates")

        let membranedata = membraneCoordinates(coordinates!, radii!, probeRad, unitcell )

        let packedcoords = membranedata.0
        let imgdata = membranedata.1
        let imgradii = membranedata.2
        
        //let packedcoords = unitcell.packcell(coordinates!)

        // atoms are red if not moved, green otherwise

        var rgb = [Vector]() 

        for (orig,tran) in zip(coordinates!,packedcoords) {
            if orig == tran {
                rgb.append(Vector([1.0,0.0,0.0]))
            }
            else {
                rgb.append(Vector([0.0,1.0,0.0]))
            }
        }

        do {
            try exportAtoms( "OUTPUT/membraneAtoms.SPH.txt", packedcoords, radii!, rgb )
        }
        catch {
            print("could not export atoms to pymol")
        }

        // apply buffer 

        /*

        let imgdata = unitcell.buildimages(packedcoords) 

        unitcell.setimages( imgdata.0, imgdata.1, imgdata.2 )

        let imgrgb = Array( repeating:Vector([1.0,1.0,0.0]), count:unitcell.imagecoords!.count) 

        var imgwithidx = [(Int,Int)]()

        for (aidx,imgindices) in unitcell.imageindices! {
            for imgidx in imgindices {
                imgwithidx.append( (imgidx,aidx ))
            }
        }

        let imginorder = imgwithidx .sorted { $0.0 < $1.0 }

        let imgradii = imginorder .map { radii![$0.1] }
        */

        let imgrgb = Array( repeating:Vector([1.0,1.0,0.0]), count:unitcell.imagecoords!.count)

        do {
            try exportAtoms( "OUTPUT/membraneImages.SPH.txt", unitcell.imagecoords!, imgradii, imgrgb )
        }
        catch {
            print("could not export atoms to pymol")
        }

        var usecoordinates = Array(coordinates!)
        var useradii = Array(radii!)

        usecoordinates = packedcoords + imgdata.0
        useradii = useradii + imgradii

        var surfdata = generateSurfaceProbes( coordinates:usecoordinates, radii:useradii, probeRadius:probeRad, 
                    levelspacing:0.5, minoverlap:0.5, numthreads:10, 
                    skipCCWContours:true, unitcell:unitcell, debugAXES:nil)

        var probes = surfdata.0 

        let procprobedata = processMembraneProbes( probes, probeRad, unitcell)

        let keepprobes = procprobedata.1
        let bufferProbes = procprobedata.2

        var str = ""

        for probe in keepprobes {
            if !probe.singleton {
                        str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 1.0 1.0\n"
                    }
                    else {
                        str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 0.0 0.0\n"
                    }

        }

        do {
                try str.write(to: URL(fileURLWithPath:"OUTPUT/probes.txt"), atomically: true, encoding: String.Encoding.utf8)
        } catch {
                print("could not write probes file !")
        }

        str = ""

        for probe in bufferProbes {
            
        str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 1.0 0.0\n"
                    

        }

        do {
                try str.write(to: URL(fileURLWithPath:"OUTPUT/imageprobes.txt"), atomically: true, encoding: String.Encoding.utf8)
        } catch {
                print("could not write probes file !")
        }

        let useprobes = procprobedata.0

        var tridata:([Vector],[Vector],[[Int]])?

        do {
            tridata = try generateTriangulation( probes:useprobes, probeRadius:probeRad, gridspacing:0.25, 
            densityDelta:0.1, densityEpsilon:0.1, isoLevel:1.0, numthreads:10, mingridchunk:20 ) 
        }
        catch {
            print("triangulation code failed !")
            exit(0)
        }

   // write out in obj format

    var VERTICES = tridata!.0 
    var NORMALS = tridata!.1
    var FACES = tridata!.2

    var url = URL(fileURLWithPath: "./membrane_init.obj")

    var outstr = ""

    for vertex in VERTICES {
        outstr += "v \(vertex.coords[0]) \(vertex.coords[1]) \(vertex.coords[2])\n"
    }

    for normal in NORMALS {
        outstr += "vn \(normal.coords[0]) \(normal.coords[1]) \(normal.coords[2])\n"
    }
    
    for face in FACES {
        outstr += "f \(face[0]+1) \(face[1]+1) \(face[2]+1)\n"
    }
    

    do {
        try outstr.write(to: url, atomically: true, encoding: String.Encoding.utf8)
    } catch {
        print("error writing file \(url)")
    }

    let procmembranetri = processMembraneTri( VERTICES, NORMALS, FACES, unitcell ) 

    VERTICES = procmembranetri.0 
    NORMALS = procmembranetri.1
    FACES = procmembranetri.2

    url = URL(fileURLWithPath: "./membrane_proc.obj")

    outstr = ""

    for vertex in VERTICES {
        outstr += "v \(vertex.coords[0]) \(vertex.coords[1]) \(vertex.coords[2])\n"
    }

    for normal in NORMALS {
        outstr += "vn \(normal.coords[0]) \(normal.coords[1]) \(normal.coords[2])\n"
    }
    
    for face in FACES {
        outstr += "f \(face[0]+1) \(face[1]+1) \(face[2]+1)\n"
    }
    

    do {
        try outstr.write(to: url, atomically: true, encoding: String.Encoding.utf8)
    } catch {
        print("error writing file \(url)")
    }






        /*
        
        let radiiVec = Vector(radii!)

        

        var allcoords = [Vector]()

        allcoords += packedcoords

        allcoords += unitcell.imagecoords!

        var allRadii = [Double]()

        allRadii += radii!

        allRadii += imgradii

        let allRadiiVec = Vector(allRadii)

        var atomcoords = [Double]() 

        for vec in allcoords {
            atomcoords += vec.coords 
        }

        let atompos = Matrix<Double>([allcoords.count,3], content:atomcoords )

        var PROBES = [Probe]()

        for axis in [AXES.X, AXES.Y, AXES.Z] {

            var axislabel = ["X", "Y", "Z"][axis.rawValue]


            var axiscoords = Vector( allcoords .map { $0.coords[axis.rawValue] } )

            let lowermin = axiscoords - allRadiiVec         
            
            let lowlim = lowermin.coords.min()! - probeRad 

            let levelspacing = 0.5

            var minaxiscoord = ceil(abs(lowlim)/levelspacing) * levelspacing 

            if lowlim < 0 {
                minaxiscoord = -minaxiscoord
            }

            //let layeridx = Int((30.0 - minaxiscoord)/levelspacing)

            // do NOT feed all coords to this function, it already adds images - maybe that is a mistake
            // as a test, feed in nil unit cell 

            let atomcircleLAYERS = atomCirclesForLayers( atompos:atompos, radii:allRadii, 
                    proberad:probeRad, minaxiscoord:minaxiscoord, layerdelta:levelspacing, axis:axis, numthreads:10, unitcell:nil )

        

            let layerBits = atomcircleLAYERS.layerBits

            let allcircles = atomcircleLAYERS.objects as! [AtomCircle]

            print("\nhave \(allcircles.count) circles total")

            for layeridx in 0..<layerBits.count {

                if layerBits[layeridx] == nil {
                    print("skipping empty layer \(layeridx)")
                    continue
                }

                let layercircles = layerBits[layeridx]!.indices() .map { allcircles[$0] }

                print("\nhave \(layercircles.count) circles in layer \(layeridx)")


                var outstr = ""

                let axisToUnit = [ Vector([1.0,0.0,0.0]),  Vector([0.0,1.0,0.0]), Vector([0.0,0.0,1.0])]

                for circle in layercircles {
                    let ident = "c\(circle.atom)"
                    let cx = circle.center.coords[0]
                    let cy = circle.center.coords[1]
                    let cz = circle.center.coords[2]
                    let rad = circle.radius 
                    let px = axisToUnit[circle.axis.rawValue].coords[0]
                    let py = axisToUnit[circle.axis.rawValue].coords[1]
                    let pz = axisToUnit[circle.axis.rawValue].coords[2]
                    var r = 1.0
                    var g = 0.0
                    var b = 0.0
                    if circle.imageatom {
                        r = 0.0
                        g = 1.0
                    }
                    outstr += "\(ident) \(cx) \(cy) \(cz) \(rad) \(px) \(py) \(pz) \(r) \(g) \(b)\n"
                }

                do {
                    try outstr.write(to: URL(fileURLWithPath:"OUTPUT/membrane.\(axislabel).CIRCLE.layer\(layeridx).txt"), atomically: true, encoding: String.Encoding.utf8)
                }
                catch {
                    print("could not write circles export for layer \(layeridx)")
                }

                
            }

            let probedata = intersectingCirclesForLayers(atomcircleLAYERS, probeRadius:probeRad, numthreads:10, skipCCWContours:true )

            let probes = probedata.1 

            var str = ""
            
            var keepprobes = [Probe]()

            for probe in probes {
                let cellc = unitcell.cellcoords(probe.center)
                var outside = false
                for ax in 0..<3 {
                    if ax == unitcell.membraneaxis.rawValue {
                        continue
                    }
                    if cellc[ax] > 1.0 || cellc[ax] < 0.0 {
                        outside = true
                    }
                }
                if outside {
                    continue
                }

                keepprobes.append(probe)

                if !probe.singleton {
                    str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 1.0 1.0\n"
                }
                else {
                    str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 0.0 0.0\n"
                }
            }

            do {
                try str.write(to: URL(fileURLWithPath:"OUTPUT/probes.\(axislabel).txt"), atomically: true, encoding: String.Encoding.utf8)
            } catch {
                print("could not write probes file !")
            }

            PROBES += keepprobes 
        }

        // add buffer for probes ; 

        let probecoords = PROBES .map { $0.center }

        let probeimgdata = unitcell.buildimages( probecoords )

        var str = ""

        for center in probeimgdata.0 {
            str += "\(center.coords[0]) \(center.coords[1]) \(center.coords[2]) \(probeRad) 1.0 1.0 0.0\n"
        }

        do {
            try str.write(to: URL(fileURLWithPath:"OUTPUT/imgprobes.txt"), atomically: true, encoding: String.Encoding.utf8)
        } catch {
            print("could not write image probes file !")
        }

        let probeinvimg = probeimgdata.2 

        var bufferPROBES = [Probe]() 

        for imgidx in probeinvimg.keys {
            let probidx = probeinvimg[imgidx]!
            let probe = PROBES[probidx]
            let bufprobe = Probe(center:probeimgdata.0[imgidx], radius:probe.proberadius, atoms:probe.atoms, 
                singleton:probe.singleton, clockwise:probe.clockwisecontour)
            bufferPROBES.append(bufprobe)


        }

        let totPROBES = PROBES + bufferPROBES

        var tridata:([Vector], [Vector], [[Int]])?

        do {
            tridata = try generateTriangulation( probes:totPROBES, probeRadius:probeRad, gridspacing:0.25, 
                densityDelta:0.1, densityEpsilon:0.1, isoLevel:1.0, numthreads:10, axis:AXES.Z)
        }
        catch {
            print("triangulation failed")
        }
    
    // write out in obj format

    let VERTICES = tridata!.0 
    let NORMALS = tridata!.1
    let FACES = tridata!.2

    let url = URL(fileURLWithPath: "./membrane.obj")

    var outstr = ""

    for vertex in VERTICES {
        outstr += "v \(vertex.coords[0]) \(vertex.coords[1]) \(vertex.coords[2])\n"
    }

    for normal in NORMALS {
        outstr += "vn \(normal.coords[0]) \(normal.coords[1]) \(normal.coords[2])\n"
    }
    
    for face in FACES {
        outstr += "f \(face[0]+1) \(face[1]+1) \(face[2]+1)\n"
    }
    

    do {
        try outstr.write(to: url, atomically: true, encoding: String.Encoding.utf8)
    } catch {
        print("error writing file \(url)")
    }

    
    // identify triangles that cross unit cell boundaries, match them?
   
    var inelems = [Int]()
    var outelems = [Int]()
    var crosselems = [Int]()

    let vertexcellcoords = VERTICES .map { unitcell.cellcoords($0) }

    for (fidx,face) in FACES.enumerated() {
        var inside = 0
        var outside = 0
        for vertex in face .map { vertexcellcoords[$0] } {
            var out = false

            for ax in 0..<3 {
                if ax == unitcell.membraneaxis.rawValue { continue }
                if vertex.coord[ax] < 0.0 || vertex.coord[ax] > 1.0 {
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

    print("\n\n# faces = \(FACES.count) , # inside = \(inelems.count) , # outside = \(outelems.count) , # crossing = \(crosselems.count)")
    
    var crossvertices = Set<Int>()

    // delete crossover elements with at 'right' or 'top' of unit cell
    // in brief, keep any element with all vertex coordinates less than 1.0

    var keepelems = Array(inelems) 

    for fidx in crosselems {
        var keep = true 
        for vertcoord in FACES[fidx] .map { vertexcellcoords[$0] }
        
        for (aidx,coord) in enumerate(vertcoord.coords) {
            if aidx == unitcell.membraneaxis.rawValue { continue }
            if coord > 1.0 {
                keep = false
                break
            }
        }

        if keep {
            keepelems.append(fidx)
        }
    }




        
        var centers = [Double]()
        var cradii = [Double]()

        for circle in layercircles {
            centers += circle.center.coords
            cradii.append(circle.radius)
        }

        let circleCenters = Matrix<Double>([layercircles.count,3], content:centers )
        let circleRadii = Matrix<Double>([layercircles.count], content:cradii )

        var pairs:[[Int]]?

        do {
            let dists = try cdist(circleCenters,circleCenters)
            try dists.setdiagonal(1.0e12)
            let sumRadii = try circleRadii.addTranspose(circleRadii)

            let mask = try Mask.compare(dists, sumRadii) { $0 < $1 }
            pairs = mask.nonzero()
        }
        catch {
            print("exception in intersectCircles for layer \(layeridx)")
            
        }


        for p in pairs! {
            if p[0] < p[1] {
                intersectAtomCircles( layercircles[p[0]], layercircles[p[1]] )

            }
        }


            // make contour for this layer

            var contours = [Contour]()

            while true {
                var cont:Contour?
                do {
                    cont = try Contour(layercircles)
                }
                catch {
                    if (error as! ContourError) == ContourError.noInitialArc {
                        print("no initial contour")
                        break
                    }
                    print("exception in Contour : \(error)")
                    break
                    
                }
                
                if cont == nil {
                    print("nil contour")
                    
                }

                contours.append(cont!)

            }

            print("\nhave \(contours.count) contours in layer \(layeridx)")

            for cont in contours {
                print("\thave \(cont.arcsInOrder.count) arcs in contour")
            }

            outstr = ""

            for (cidx,cont) in contours.enumerated() {
                for (arcidx,arc) in cont.arcsInOrder.enumerated() {
                    let ctr = arc.parentcircle.center
                    let pstart = arc.pstart
                    let pend = arc.pend 
                    let psx = pstart.coords[0]
                    let psy = pstart.coords[1]
                    let psz = pstart.coords[2]
                    let pex = pend.coords[0]
                    let pey = pend.coords[1]
                    let pez = pend.coords[2]

                    
                    let cx = ctr.coords[0]
                    let cy = ctr.coords[1]
                    let cz = ctr.coords[2]
                    let rad = arc.parentcircle.radius 
                    let px = axisToUnit[arc.parentcircle.axis.rawValue].coords[0]
                    let py = axisToUnit[arc.parentcircle.axis.rawValue].coords[1]
                    let pz = axisToUnit[arc.parentcircle.axis.rawValue].coords[2]

                    let ident = "c_\(cidx)_\(arcidx)"

                    var r = 0.0
                    var g = 0.0
                    var b = 1.0

                    if arcidx == 0 {
                        g = 1.0
                        b = 0.0
                    }
                    else if arcidx == cont.arcsInOrder.count - 1 {
                        r = 1.0
                        b = 0.0
                    }

                    outstr += "\(ident) \(cx) \(cy) \(cz) \(psx) \(psy) \(psz) \(pex) \(pey) \(pez) \(px) \(py) \(pz) \(r) \(g) \(b)\n"

                }
            }

            do {
                try outstr.write(to: URL(fileURLWithPath:"membrane.CONTOURS.layer\(layeridx).txt"), atomically: true, encoding: String.Encoding.utf8)
            }
            catch {
                print("could not write contours export for layer \(layeridx)")
            }
 
        }
    



       
        // OK, at this point do all the layers + probes

        let probedata = intersectingCirclesForLayers(atomcircleLAYERS, probeRadius:probeRad, numthreads:10, skipCCWContours:true )
        
        
        let probes = probedata.1 

        var str = ""
        

        for probe in probes {
            if !probe.singleton {
                str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 1.0 1.0\n"
            }
            else {
                str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 0.0 0.0\n"
            }
        }

        do {
            try str.write(to: URL(fileURLWithPath:"probes.txt"), atomically: true, encoding: String.Encoding.utf8)
        } catch {
            print("could not write probes file !")
        }

        
       */ 

    }

    /*
    func testCirclesArcs() {

        // make some circles with simple geometry, test intersections

        let PI = acos(-1.0)

        var axis = AXES.X 

        let centerA = Vector([0.0,0.0,0.0])
        let centerB = Vector([0.0,2.0,0.0])
        let centerC = Vector([0.0,1.0,2.0*sin(60.0*(PI/180.0))])

        var circleA = AtomCircle(0, centerA, 1.25,  axis )
        var circleB = AtomCircle(1, centerB, 1.25,  axis )
        var circleC = AtomCircle(2, centerC, 1.25,  axis )

        print("Before intersection : ")
        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")

        // intersect A and B

        intersectAtomCircles( circleA, circleB )

        print("\nAfter intersecting A and B : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")

        // intersect C with A and B

        intersectAtomCircles( circleA, circleC )
        intersectAtomCircles( circleB, circleC )

        print("\nAfter second intersection, A and B with C  : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")

        // check whether vectors inside arcs

        let uAC = centerC.sub(centerA).unit()!
        var inside = circleA.exposure[0].contains(uAC)
        print("\nuAC inside A exposure = \(inside)")
        inside = circleC.exposure[0].contains(uAC)
        print("\nuAC inside C exposure = \(inside)")
        inside = circleB.exposure[0].contains(uAC)
        print("\nuAC inside B exposure = \(inside)")

        // add small circle D at upper intersection of A and B, center =  0.0, 1.0, 0.75
        // should not affect exposures

        let centerD = Vector([0.0, 1.0, 0.75])
        var circleD = AtomCircle(3, centerD, 0.25,  axis )

        intersectAtomCircles( circleA, circleD )
        intersectAtomCircles( circleB, circleD )
        intersectAtomCircles( circleC, circleD )
        
        print("\nAfter adding circle D, intersecting with existing A,B,C  : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")
        print("D: \(circleD.str())")

        
        let centerE = Vector([0.0, 1.0, -0.75])
        var circleE = AtomCircle(4, centerE, 0.25,  axis )

        intersectAtomCircles( circleA, circleE )
        print("After intersecting A and E :")
        print("A: \(circleA.str())")
        print("E: \(circleE.str())")

        intersectAtomCircles( circleB, circleE )
        print("After intersecting B and E :")
        print("B: \(circleB.str())")
        print("E: \(circleE.str())")

        //
        intersectAtomCircles( circleC, circleE )
        print("After intersecting C and E :")
        print("C: \(circleC.str())")
        print("E: \(circleE.str())")
        intersectAtomCircles( circleD, circleE )
        
        print("\nAfter adding circle E, intersecting with existing A,B,C,D  : ")

        print("A: \(circleA.str())")
        print("B: \(circleB.str())")
        print("C: \(circleC.str())")
        print("D: \(circleD.str())")    
        print("E: \(circleE.str())")

        var circles = [circleA,circleB,circleC,circleE]

        var contour:Contour?

        do {
            contour = try Contour(circles)
        }
        catch {
            print("exception in Contour")
        }

        print(contour!.str(verbose:true))

        // Regenerate circles, try pairwise

        print("pairwise intersection")

        circleA = AtomCircle(0, centerA, 1.25,  axis )
        circleB = AtomCircle(1, centerB, 1.25,  axis )
        circleC = AtomCircle(2, centerC, 1.25,  axis )
        circleD = AtomCircle(3, centerD, 0.25,  axis )
        circleE = AtomCircle(4, centerE, 0.25,  axis )

        circles = [circleA,circleB,circleC,circleD,circleE]

        for i in 0..<(circles.count-1) {
            for j in (i+1)..<circles.count {
                intersectAtomCircles(circles[i],circles[j])
            }
        }

        do {
            contour = try Contour(circles)
        }
        catch {
            print("exception in Contour")
        }

        print(contour!.str(verbose:true))

        print("pairwise intersection, change order")

        circleA = AtomCircle(0, centerA, 1.25,  axis )
        circleB = AtomCircle(1, centerB, 1.25,  axis )
        circleC = AtomCircle(2, centerC, 1.25,  axis )
        circleD = AtomCircle(3, centerD, 0.25,  axis )
        circleE = AtomCircle(4, centerE, 0.25,  axis )

        circles = [circleC,circleE,circleA,circleD,circleB]

        for i in 0..<(circles.count-1) {
            for j in (i+1)..<circles.count {
                intersectAtomCircles(circles[i],circles[j])
            }
        }

        do {
            contour = try Contour(circles)
        }
        catch {
            print("exception in Contour")
        }

        print(contour!.str(verbose:true))

    }

    
    
    
    func testLayers() throws {

        


        // Make collection of atom circles

        
        // read circle definitions from data directory
        
        var importdata:([Vector]?,[Double]?,Double?,[Double]?,[[AtomCircle]]?,[[[(Int,Vector?,Int?,Vector?,Int?)]]]?)?
        
        do {
            importdata = try importData("Tests/DATA/data.txt")
        }
        catch {
            print("exception calling importData")
            throw tensurError.testImportError
        }

        let data = importdata!

        let coordinates = data.0!
        let radii = data.1!
        let probeRad = data.2!
        let levels = data.3!
        let circlesForLevels = data.4!
        let contoursForLevels = data.5!

        var atomcoords = [Double]() 

        for vec in coordinates {
            atomcoords += vec.coords 
        }

        let atompos = Matrix<Double>([coordinates.count,3], content:atomcoords )

        let maxrad = radii.max()!

        let xcoords = coordinates.map { $0.coords[0] }
        let minX = xcoords.min()!

        let mincoord = levels.min()!

        print("min layer coord this data = \(mincoord)")
        print("min X-coord = \(minX), max rad = \(maxrad), min X-coord - maxrad - probeRad = \(minX - maxrad - probeRad)")

        let circleLayers = atomCirclesForLayers( atompos:atompos, radii:radii, 
            proberad:probeRad, minaxiscoord:mincoord, layerdelta:0.5, axis:AXES.X, numthreads:1 )
            
        // Get all atom circles with X 
        
        var atomcircles = circleLayers.edgeToFaces[edge]ects as! [AtomCircle]

        let X = 26.5

        var layercircles = atomcircles.filter { $0.center.coords[0] == X }
        
        print("have \(layercircles.count) circles for layer at \(X)")


        //print("layer circle data ---")

        //for circle in layercircles {
        //    print(circle.str())
        //}

        for i in 0..<(layercircles.count-1) {
            for j in (i+1)..<layercircles.count {
                //let d = layercircles[i].center.sub(layercircles[j].center).length()
                //if d > layercircles[i].radius + layercircles[j].radius {
                //    continue
                //}
                
                    //print("\n\nintersect \(i) and \(j)\n")
                    //print("prior intersect, arcs \(i)---")
                    //for arc in layercircles[i].exposure {
                    //    print(arc.str())
                    //}
                    //print("prior intersect, arcs \(j)---")
                    //for arc in layercircles[j].exposure {
                    //    print(arc.str())
                    //}
                
                
                intersectAtomCircles(layercircles[i],layercircles[j])
                
                    //print("after intersect, arcs \(i)---")
                    //for arc in layercircles[i].exposure {
                    //    print(arc.str())
                    //}
                    //print("after intersect, arcs \(j)---")
                    //for arc in layercircles[j].exposure {
                    //    print(arc.str())
                    //}
                
            }
        }

        // look for singleton circle

        let UP = axisUP[AXES.X.rawValue]
        let RIGHT = axisRIGHT[AXES.X.rawValue]

        //print("singleton circle")

        for circle in layercircles {
            if circle.center.coords[UP] > 5.0 && circle.center.coords[RIGHT] < 18.0 {
                print(circle.str())
                //print("exposure: \(circle.exposure.count)")
            }
        }

        var contour:Contour?

        do {
            contour = try Contour(layercircles)
        }
        catch {
            print("exception in Contour first round")
        }

        if contour != nil {
            //print("CONTOUR, first round ----- ")
            //print(contour!.str(verbose:true))
        }
        else {
            print("contour failed first round")
        }

        //print("\ncontour python ...\n")

        //let text = printPython( layercircles, contour!.arcsInOrder)

        //print(text)

        // check for removed arcs - OK working now

        
        //print("removed circles/arcs ---")

        for circle in layercircles {
            if circle.removed {
                //print("removed circle at atom \(circle.atom)")
                continue
            }
            else {
                //print("retained circle at atom \(circle.atom)")
            }
            for arc in circle.exposure {
                if arc.removed {
                    //print("removed arc")
                    //print(arc.str())
                }
            }
        }
       

        // print all arcs 

        //print("all arcs :")

        //for circle in layercircles {
        //    if circle.removed {
        //        continue 
        //    }
        //    for arc in circle.exposure {
        //        print(arc.str())
        //    }
        //}
        
        

        //var axis:AXES

        //var mincoord:Double
        //var maxcoord:Double

        //var delta:Double

        //var objects:[Any]

        //var layers:[Int]

        //var layerBits:[objectBits?]

        print("layers min, max coord = \(circleLayers.mincoord), \(circleLayers.maxcoord)")
        print("layers count = \(circleLayers.layerBits.count)")
        print("have \(circleLayers.objects.count) atom circles")
        print("\n\n")

        
 
        // Get all atomcircles inside another and remove 
    

        atomcircles = circleLayers.objects as! [AtomCircle]
        var centers = [Double]() 
        var centerradii = atomcircles.map { $0.radius }

        for atomcircle in atomcircles {
            centers += atomcircle.center.coords
        }

        print("have \(atomcircles.count) atom circles, length centers array = \(centers.count)")

        for lidx in 0..<circleLayers.layerBits.count {
            let lb = circleLayers.layerBits[lidx]
            if lb == nil {
                //print("layer \(lidx) : nil ")
            }
            else {
                let indices = lb!.indices()
                let first = indices[0]
                let x = atomcircles[first].center.coords[0]
                let y = atomcircles[first].center.coords[1]
                let z = atomcircles[first].center.coords[2]
                //print("layer \(lidx) : \(lb!.indices().count) atom circles, center at \(x), \(y), \(z) ")
            }
            
        }

        
        let surfdata = generateSurfaceProbes( coordinates:coordinates, radii:radii, probeRadius:probeRad, levelspacing:0.5, minoverlap:0.5, numthreads:10,
                         skipCCWContours:true )

        let probes = surfdata.0 

        var str = ""
        

        for probe in probes {
            if !probe.singleton {
                str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 1.0 1.0\n"
            }
            else {
                str += "\(probe.center.coords[0]) \(probe.center.coords[1]) \(probe.center.coords[2]) \(probeRad) 1.0 0.0 0.0\n"
            }
        }

        do {
            try str.write(to: URL(fileURLWithPath:"probes.txt"), atomically: true, encoding: String.Encoding.utf8)
        } catch {
            print("could not write probes file !")
        }
        
        let time1 = Date().timeIntervalSince1970

        var tridata:([Vector],[Vector],[[Int]])?

        // set axis for testing 

        let testaxis = AXES.Z

        do {
            tridata = try generateTriangulation( probes:probes, probeRadius:probeRad, gridspacing:0.4, 
            densityDelta:0.1, densityEpsilon:0.1, isoLevel:1.0, numthreads:10, axis:testaxis) 
        }
        catch {
            print("triangulation code failed !")
        }

    // write out in obj format

    let VERTICES = tridata!.0 
    let NORMALS = tridata!.1
    let FACES = tridata!.2

    let url = URL(fileURLWithPath: "./isosurface.axis.\(testaxis.rawValue).obj")

    var outstr = ""

    for vertex in VERTICES {
        outstr += "v \(vertex.coords[0]) \(vertex.coords[1]) \(vertex.coords[2])\n"
    }

    for normal in NORMALS {
        outstr += "v \(normal.coords[0]) \(normal.coords[1]) \(normal.coords[2])\n"
    }
    
    for face in FACES {
        outstr += "f \(face[0]+1) \(face[1]+1) \(face[2]+1)\n"
    }
    

    do {
        try outstr.write(to: url, atomically: true, encoding: String.Encoding.utf8)
    } catch {
        print("error writing file \(url)")
    }

    print("wrote file to path \(url)")

    let time2 = Date().timeIntervalSince1970

    print("\ntime for marching cubes = \(time2 - time1)")


        
    }
    
   
    */

    

 
}

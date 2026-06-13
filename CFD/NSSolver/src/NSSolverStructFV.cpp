#include "NSSolverStructFV.h"
#include "Flux_Inviscid.h"
#include "Constants.h"
#include "TK_Exit.h"
#include "GlobalDataBase.h"

namespace PHSPACE
{
NSSolverStructFV::NSSolverStructFV() {}
NSSolverStructFV::~NSSolverStructFV() {}

void NSSolverStructFV::InviscidFlux(Grid* gridIn) {
    const std::string scheme_name = GlobalDataBase::GetStrParaFromDB("inviscidScheme");
    
    const int scheme_id = GetSchemeID(scheme_name);
    
    switch (scheme_id) {
        case ISCHEME_ROE:
        case ISCHEME_HLLE:
            //其他格式暂不实现，保持框架结构完整
            break;
        case ISCHEME_GKS:
            GKSFluxStruct(gridIn);
            break;
        default:
            TK_Exit::ExceptionExit("NSSolverStructFV::InviscidFlux: 未知通量格式 " + scheme_name);
    }
}
//GKS主函数：框架级接口
void NSSolverStructFV::GKSFluxStruct(Grid* gridIn) {
    if (gridIn == nullptr) {
        TK_Exit::ExceptionExit("GKSFluxStruct:输入网格为空");
    }

    std::cout << "[GKS框架]通量计算链路激活，网格地址: " << gridIn << std::endl;
    //预留数据流向框架（后续扩展区）
}

void NSSolverStructFV::ComputeGKSFluxCore(
    RDouble rhoL, RDouble uL, RDouble vL, RDouble wL, RDouble eL,
    RDouble rhoR, RDouble uR, RDouble vR, RDouble wR, RDouble eR,
    RDouble nx, RDouble ny, RDouble nz,
    RDouble* flux
) {
    if (flux != nullptr) {
        for (int i = 0; i < 5; ++i) {
            flux[i] = 0.0;//临时初始化，后替换为实际公式
        }
    }
}
}
# Claude에서 이어서 수정하기

> 2026-10-14 갱신: 이 인계 이후 10월 13일 리뷰 지적사항을 반영했습니다. 현재 상태는
> `internal/REVIEW_codex_package_20261011.md`의 Resolution 절과 `CLAUDE.md`를 기준으로 보십시오.
> 아래 내용은 인계 시점의 기록입니다.

이번 커밋은 저자가 요청한 로컬 최종 수정본의 GitHub 인계입니다. 논문 수치나
문장을 추가로 고치지 않고 2026-10-03 분석 / 2026-10-11 전달본을 반영했습니다.
원격에 먼저 올라온 10월 13일 리뷰는 그대로 보존했습니다.

## 먼저 읽을 파일
1. `MANUSCRIPT_TED_TRAP_v5_MASTER.md`: 이번에 반영한 최신 원고, 유일한 MD master.
2. `submission/candidate_20261011_review/`: DOCX/PDF, 본문 및 보충 표·그림, XLSX.
3. `internal/REVIEW_codex_package_20261011.md`: 아직 해결하지 않은 Claude 리뷰.
4. `submission/candidate_20261011_review/FINAL_RESULTS_INSIGHTS_KO.md`: 결과 설명.
5. `submission/candidate_20261011_review/reproducibility/README.md`: 계산 코드 안내.

## 9월 원격 버전에서 달라진 내용
- FinnGen GD와 strict Graves ophthalmopathy 추가: 다섯 질병 자료, 서로 겹치는
  FinnGen 세 endpoint는 독립 코호트가 아님.
- 관련 유전자 비교, 24개 유전자 지역 공위치화와 민감도 결과, TSHR 대립유전자 점검.
- TSHR/IGF1R 5개 조직 × 5개 질병의 50개 계획 비교, 실제 계산 20개 / 불가 30개.
- Table 3A–C 통합, Tables S1–S3 및 Figures S1–S2 유지, 상세데이터 XLSX 하나.
- 원고·표·그림·보충자료 통합본과 원본 코드/집계 결과를 함께 인계.

## 다음 수정의 출발점
10월 13일 리뷰의 초록/결론 표현, 공위치화 불일치 노출, 지방조직 결과의 독립성,
FinnGen GD IGF1R 추정치 해석, 그림 범례/표 서식/보고사항, 원자료 기반 재현성
검토를 이어가십시오. 이 인계는 해당 지적을 해결했다고 주장하지 않습니다.
리뷰의 각 항목은 실제 파일과 데이터로 확인한 다음 반영해야 합니다.
저자는 Claude에서 추가 수정할 예정이므로 이번 커밋에는 새 분석이나 문장 수정이 없습니다.

## 파일 보존 및 검증
- 9월 원격 candidate는 `archives/submission_candidates/candidate_20260911_review/`로 이동.
- 이전 master/README/CLAUDE/audit/checklist는 `internal/handoff_before_20261014/`에 보존.
- 원본 데이터 추출물, 개인별 참조 유전형, 실행 로그는 로컬에 유지하고 커밋에서 제외.
  원 논문용 DOCX/PDF/PNG/XLSX는 바이트 동일; 문서 및 코드의 로컬 사용자 경로만 익명화.
- 파일 해시와 조직 MR 식 검사는 `python scripts/audit_paper1_integrity.py`.
  이 검사는 원격에 없는 영역 통계의 재계산이나 리뷰 완료를 대신하지 않습니다.
- 원본 실행 코드는 workspace 상대경로를 사용합니다. 원자료 및 경로 복원이 필요하며,
  기존 결과 CSV를 감사에 맞추어 수정하지 마십시오.
- 최종 학술지/저자 승인/선언 문구/cover letter는 후속 검토에서 확정하십시오.
